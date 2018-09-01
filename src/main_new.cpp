#include <bits/stdc++.h>
#include <random>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>

using namespace std;

default_random_engine random_engine;
default_random_engine random_engine2;
const double eps = 1e-9;
const bool isOutputMode = 1;
const int seed = 23335;
const int inf = 1 << 29;
const int REGION_NUM = 30;
const int DRIVER_NUM = 105;
const int REQUEST_NUM = 400000;
const int TIME_STEP = 60;
const int TOTAL_TIME_INTERVAL = 24 * (3600 / TIME_STEP);
const string mode = "file"; 			// "file" or "synthetic"

const int TOTAL_TIMES = 50;
const int TIMES_FOR_ONE_DISTRIBUTION = 20;

const string Valuemode = "no"; 			// "random" or "no"
const string Rejectionmode = "no"; 		// "random" or "strong" or "no"
const string KMmode = "without future"; // "with future" or " "without future""
// const string stage1mode = "first-fit"; 	// "first-fit" or "heuristic"
// const string stage2mode = "KM"; 		// "normal" or "MDP" or "KM"
const string stage1mode = "heuristic"; 	// "first-fit" or "heuristic"
const string stage2mode = "MDP"; 		// "normal" or "MDP" or "KM"
const int BIG_DRIVER_NUM = 45000;

// manual settings for "synthetic" mode
const int DEBUG = 1; // 0 means no
const int MIN_DIST = 100;
const int MAX_DIST = 1000;
const int MIN_TRIP_TIME = 5 * 60;
const int MAX_TRIP_TIME = 40 * 60;
//const int MIN_WAITING_TIME = 5 * 60;
//const int MAX_WAITING_TIME = 30 * 60;
int driver_num = 100;
int total_driver_num = 100;
int regions = 21;
int pre_num = 100000;
int ond_num = 200000;

// manual settings for "file" mode
const int REAL_DRIVER_NUM = 50;
const double default_ratio = 10;
const double pre_ratio = 0.05;
const double ond_ratio = 1.00;
double driver_ratio = 1.;
int days = 10;

// heuristic settings
const double a = 1.;
const double b = 1.;
const double k = 1.;
uniform_real_distribution<double> U(0.0, 1.0);

double VALUE1 = 0;
int COUNT = 0;
int times;

double calc_heuristic_score_function(double delta) {
	double alpha = U(random_engine2);
	double rank = exp(k * alpha);
    VALUE1 += delta;
    COUNT ++;
	return a * rank + b * delta;
}

struct SimpleRequest {
    int order, start_time, trip_time, waiting_time, advanced_time;
    int source, destination, isPre;

    SimpleRequest() {}

    SimpleRequest(const SimpleRequest &a) {
        order = a.order;
        start_time = a.start_time;
        trip_time = a.trip_time;
        waiting_time = a.waiting_time;
        advanced_time = a.advanced_time;
        source = a.source;
        destination = a.destination;
        isPre = a.isPre;
    }

    void print() {
        printf("order:%d\tstart time:%d\ttrip time:%d\tadvanced time:%dis Pre%d\n", order, start_time, trip_time, advanced_time, isPre);
    }
};

int dist[REGION_NUM][REGION_NUM];
int driver_time[BIG_DRIVER_NUM], driver_src[BIG_DRIVER_NUM];
int satisfied_ond = 0, satisfied_pre = 0;
int Count[REGION_NUM][REGION_NUM][TOTAL_TIME_INTERVAL * 2];
int CountValue[REGION_NUM][REGION_NUM][TOTAL_TIME_INTERVAL * 2];
double Prob[REGION_NUM][REGION_NUM][TOTAL_TIME_INTERVAL * 2][DRIVER_NUM];
double initProb[REGION_NUM][REGION_NUM][TOTAL_TIME_INTERVAL * 2][DRIVER_NUM];
double V[REGION_NUM][REGION_NUM][TOTAL_TIME_INTERVAL * 2];

struct ActiveDriver {
    int driver, time, number;
    int u, v, length, finished_work;
    int type; // 0 means ilde, 1 means full
    vector<SimpleRequest> missions;
    ActiveDriver(int driver = -1) : driver(driver), time(driver_time[driver]), number(driver), u(driver_src[driver]), v(-1), length(-1), finished_work(0), type(0), missions() {}
} drivers[DRIVER_NUM];

struct Event {
    int type; // 0 means an on-demnad order; 1 means a driver;
    int time, number; // number means driver-number / request-number

    Event() {
        type = -1;
        time = -1;
        number = -1;
    }

    Event(const Event &e) {
        type = e.type;
        time = e.time;
        number = e.number;
    }
};

inline bool operator<(const Event &e1, const Event &e2) {
    return e1.time > e2.time || (e1.time == e2.time && e1.type > e2.type);
}

void printMissions(int driver) {
    printf("---Driver %d has %d missions. Finished %d missions.\n", driver, int(drivers[driver].missions.size()), drivers[driver].finished_work);
    printf("Driver from %d to %d, length: %d\n", drivers[driver].u, drivers[driver].v, drivers[driver].length);
    for (int j = 0; j < drivers[driver].missions.size(); j ++) {
        SimpleRequest sr = drivers[driver].missions[j];
        printf("Mission %d\tRequest %02d\tsource: %d\tdestination: %d\tstart time: %d\toccur time: %d\ttrip time: %d\n", j + 1, sr.order, sr.source, sr.destination, sr.start_time, max(0, sr.start_time - sr.advanced_time), sr.trip_time);
        if (j > 0) printf("dist:%d\n", dist[drivers[driver].missions[j - 1].destination][sr.source]);
    }
    puts("");
}

void printRequest(const SimpleRequest &sr) {
    printf("Request %02d\tsource: %d\tdestination: %d\tstart time: %d\toccur time: %d\ttrip time: %d\n", sr.order, sr.source, sr.destination, sr.start_time, max(0, sr.start_time - sr.advanced_time), sr.trip_time);
    puts("");
}

vector<SimpleRequest> pre_requests, ond_requests, offline_requests;
vector<SimpleRequest> pre_requests_days[20], ond_requests_days[20];
int assigned_driver[REQUEST_NUM];

int calc_time(const int &l1, const int &l2) {
    return dist[l1][l2];
}

inline int remainedMission(int driver) {
    return int(drivers[driver].missions.size()) - drivers[driver].finished_work;
}

int total_pre_value = 0, total_ond_value = 0;
int total_satisfied_pre_value = 0, total_satisfied_ond_value = 0;
bool isChosenDriver[BIG_DRIVER_NUM];

double g[REGION_NUM][REGION_NUM][TOTAL_TIME_INTERVAL];
map<pair<int, int>, pair<int, double>> Map;
double pre_computed_table[REGION_NUM][TOTAL_TIME_INTERVAL][REGION_NUM][TOTAL_TIME_INTERVAL];
int time_stamp = 0;

inline bool cmp_score(const pair<int, pair<double, double>> &a, const pair<int, pair<double, double>> &b) {
    return a.second.first + a.second.second > b.second.first + b.second.second;
}

double calc_score_value(const int &v, const int &t, const int &v2, const int &t2) {
//	if (check) cout << v << ' ' << t << endl;
    if (v == v2 && t == t2) return 0;
    if (Map[make_pair(v, t)].first == time_stamp) {
        return Map[make_pair(v, t)].second;
    }
    vector<pair<int, pair<double, double>>> available_regions;
    for (int d = 0; d < regions; d ++) {
        if ((d == v2 ? 0 : dist[d][v2]) <= t2 - t - dist[v][d]) {
            available_regions.push_back(make_pair(d, make_pair(V[v][d][t], calc_score_value(d, t + dist[v][d], v2, t2))));
        }
    }
    sort(available_regions.begin(), available_regions.end(), cmp_score);
    double d_star = 0;
    for (auto tmp : available_regions) {
        d_star = max(d_star, tmp.second.second);
    }
    // if (t == 0)
    //     for (auto tmp : available_regions)
    //         cerr << tmp.second.first << ' ' << tmp.second.second << endl;
    if (v == v2 && t + 1 <= t2) {
        d_star = max(d_star, calc_score_value(v, t + 1, v2, t2));
    } else if (v != v2 && dist[v][v2] <= t2 - t - 1) {
        d_star = max(d_star, calc_score_value(v, t + 1, v2, t2));
    }
    // if (t == 0) cerr << d_star << endl;
    int j = -1;
    for (int i = 0; i < available_regions.size(); i ++) {
        if (available_regions[i].second.first + available_regions[i].second.second > d_star) {
            j = i;
        }
    }
    double p = 1.0, F = 0.0;
    for (int i = 0; i <= j; i ++) {
        int ai = available_regions[i].first;
        F += p * Prob[v][ai][t][1] * (V[v][ai][t] + calc_score_value(ai, t + dist[v][ai], v2, t2));
        p *= (1.0 - Prob[v][ai][t][1]);
		if (p > 1. + eps) {
			// cout << v << ' ' << ai<< ' ' << t << endl;
			// cout << Prob[v][ai][t][1] << endl;
			assert(false);
		}
    }
	if (F < -eps) {
	    // for (int i = 0; i <= j; i ++) {
	    //     int ai = available_regions[i].first;
		// 	cout << "prob " << p << " a[i] " << ai << endl;
		// 	cout << "V " << V[v][ai][t] << " f() " << calc_score_value(ai, t + dist[v][ai], request) << endl;
	    //     F += p * Prob[v][ai][t][1] * (V[v][ai][t] + calc_score_value(ai, t + dist[v][ai], request));
	    //     p *= 1.0 - Prob[v][ai][t][1];
	    // }
		assert(false);
	}
    F += p * d_star;
    Map[make_pair(v, t)].first = time_stamp;
    Map[make_pair(v, t)].second = F;
    return F;
}

void load_data(const string &mode) {
    if (mode == "file") {
        freopen("../data/distance.txt", "r", stdin);
        regions = 21;
        for (int i = 0; i < regions; i ++)
            for (int j = 0; j < regions; j ++) {
                scanf("%d", &dist[i][j]);
                if (dist[i][j] <= 0) {
                    dist[i][j] = random_engine() % (MAX_DIST - MIN_DIST) + MIN_DIST;
                }
                dist[i][j] /= TIME_STEP;
            }
        freopen("../data/driver.txt", "r", stdin);
        scanf("%d", &total_driver_num);
		driver_num = total_driver_num;
        for (int i = 0; i < driver_num; i ++) {
            scanf("%d%d", &driver_time[i], &driver_src[i]);
			driver_time[i] /= TIME_STEP;
        }
        freopen("../data/all_the_requests.txt", "r", stdin);
        pre_num = 0;
        ond_num = 0;
        int total_num;
        scanf("%d", &total_num);
        for (int i = 0; i < total_num; i ++) {
            SimpleRequest sr;
            scanf("%d%d%d%d%d", &sr.start_time, &sr.trip_time, &sr.waiting_time, &sr.source, &sr.destination);
            if (sr.trip_time <= 180 || sr.start_time / TIME_STEP <= 0) continue;
            sr.start_time /= TIME_STEP;
            sr.trip_time /= TIME_STEP;
//            sr.trip_time = dist[sr.source][sr.destination];
            sr.waiting_time /= TIME_STEP;
            sr.waiting_time = 0;
            long double random_value = 1.0 * (random_engine() % 10000) / 10000;
            int day = random_engine() % days;
			ond_requests_days[day].push_back(sr);
			offline_requests.push_back(sr);
        }
    } else if (mode == "synthetic") {
        for (int i = 0; i < regions; i ++)
            for (int j = 0; j < regions; j ++) {
                if (i <= j) {
                    dist[i][j] = random_engine() % (MAX_DIST - MIN_DIST) + MIN_DIST;
                } else {
                    dist[i][j] = dist[j][i];
                }
                dist[i][j] /= TIME_STEP;
            }
        for (int i = 0; i < driver_num; i ++) {
            driver_src[i] = random_engine() % regions;
            drivers[i] = ActiveDriver(i);
        }
        for (int i = 0; i < pre_num; i ++) {
            int day = random_engine() % days;
            SimpleRequest sr;
            sr.order = i;
            sr.start_time = random_engine() % TOTAL_TIME_INTERVAL; // 00:00 - 23:59
            sr.trip_time = (random_engine() % (MAX_TRIP_TIME - MIN_TRIP_TIME) + MIN_TRIP_TIME) / TIME_STEP; // 5min - 45min
            sr.source = random_engine() % regions;
            sr.destination = random_engine() % regions;
            sr.isPre = 1;
            pre_requests_days[day].push_back(sr);
            if (day == days - 1) {
                pre_requests.push_back(sr);
            }
        }
        for (int i = 0; i < ond_num; i ++) {
            int day = random_engine() % days;
            SimpleRequest sr;
            sr.order = i;
            sr.start_time = random_engine() % TOTAL_TIME_INTERVAL; // 00:00 - 23:59
            sr.trip_time = (random_engine() % (MAX_TRIP_TIME - MIN_TRIP_TIME) + MIN_TRIP_TIME) / TIME_STEP; // 5min - 45min
            // sr.waiting_time = random_engine() % (MAX_WAITING_TIME - MIN_WAITING_TIME)+ MIN_WAITING_TIME; // 5min - 30min
            sr.source = random_engine() % regions;
            sr.destination = random_engine() % regions;
            sr.isPre = 0;
            ond_requests_days[day].push_back(sr);
            if (day == days - 1) {
                ond_requests.push_back(sr);
            } else {
                offline_requests.push_back(sr);
            }
        }
    }
    for (int day = 0; day < days; day ++) {
        for (int i = 0; i < ond_requests_days[day].size(); i++) {
            SimpleRequest sr = ond_requests_days[day][i];
            Count[sr.source][sr.destination][sr.start_time] ++;
        }
        for (int i = 0; i < regions; i ++)
            for (int j = 0; j < regions; j ++)
                for (int k = 0; k < TOTAL_TIME_INTERVAL; k ++)
                    if (Count[i][j][k] >= 1) {
                        for (int l = 1; l <= Count[i][j][k]; l++) {
                            Prob[i][j][k][l] += 1.0 / days / default_ratio;
                            initProb[i][j][k][l] += 1.0 / days / default_ratio;
							assert(Prob[i][j][k][l] <= 1. + eps && Prob[i][j][k][l] >= 0);
						}
                        Count[i][j][k] = 0;
                    }
    }
    for (int day = 0; day < days; day ++) {
        for (int i = 0; i < ond_requests_days[day].size(); i++) {
            SimpleRequest sr = ond_requests_days[day][i];
            Count[sr.source][sr.destination][sr.start_time]++;
            CountValue[sr.source][sr.destination][sr.start_time] += sr.trip_time;
        }
    }
    for (int i = 0; i < regions; i ++)
        for (int j = 0; j < regions; j ++)
            for (int k = 0; k < TOTAL_TIME_INTERVAL; k ++)
                if (Count[i][j][k] >= 1) {
                    V[i][j][k] += 1.0 * CountValue[i][j][k] / Count[i][j][k];
                }
    for (int i = 0; i < ond_requests.size(); i ++) {
        assigned_driver[i] = -1;
    }
    printf("We have %d drivers.\n", driver_num);
    printf("We have %d offline requests.\n", offline_requests.size());

    time_t start_t, end_t;
    time(&start_t);
	memset(pre_computed_table, 0, sizeof(pre_computed_table));
    for (int v2 = 0; v2 < regions; v2 ++) {
		cout << v2 << endl;
        for (int t2 = 0; t2 < TOTAL_TIME_INTERVAL; t2 ++) {
            time_stamp ++;
			Map.clear();
            for (int v = 0; v < regions; v ++)
                for (int t = 0; t < t2; t ++) {
                    pre_computed_table[v][t][v2][t2] = calc_score_value(v, t, v2, t2);
                }
        }
	}
    time(&end_t);
    Map.clear();
	cout << start_t << ' ' << end_t << endl;
    printf("After %d seconds, we have computed all the values in Lool-Up Table.\n", end_t - start_t);
}

struct EdgeList {
    int size;
    int last[10000000];
    int succ[10000000], other[10000000], flow[10000000], cost[10000000];
    void clear(int n) {
        size = 0;
        std::fill(last, last + n, -1);
    }
    void add(int x, int y, int c, int w) {
        succ[size] = last[x];
        last[x] = size;
        other[size] = y;
        flow[size] = c;
        cost[size++] = w;
    }
} e;

int n, source, target;
int previous[10000000];

void add_edge(int x, int y, int c, int w) {
	// cout << x << ' ' << y << ' ' << c << ' ' << w << endl;
    e.add(x, y, c, w);
    e.add(y, x, 0, -w);
}

bool augment() {
    static int dist[10000000], occur[10000000];
    std::vector<int> queue;
    std::fill(dist, dist + n, -1);
    std::fill(occur, occur + n, 0);
    dist[source] = 0;
    occur[source] = true;
    queue.push_back(source);
    for (int head = 0; head < (int)queue.size(); ++head) {
        int x = queue[head];
        for (int i = e.last[x]; ~i; i = e.succ[i]) {
            int y = e.other[i];
            if (e.flow[i] && dist[y] < dist[x] + e.cost[i]) {
                dist[y] = dist[x] + e.cost[i];
                previous[y] = i;
                if (!occur[y]) {
                    occur[y] = true;
                    queue.push_back(y);
                }
            }
        }
        occur[x] = false;
    }
    return dist[target] > 0;
}

int max_flow() {
    std::pair<int, int> answer = std::make_pair(0, 0);
    while (augment()) {
        int number = INT_MAX;
        for (int i = target; i != source; i = e.other[previous[i] ^ 1]) {
            number = std::min(number, e.flow[previous[i]]);
        }
        answer.first += number;
        for (int i = target; i != source; i = e.other[previous[i] ^ 1]) {
            e.flow[previous[i]] -= number;
            e.flow[previous[i] ^ 1] += number;
            answer.second += number * e.cost[previous[i]];
        }
    }
    return answer.second;
}

double calc_offline() {
	n = pre_requests.size() * 2 + driver_num * 2 + 2;
	source = n - 1;
	target = n - 2;
	e.clear(n);
	for (int i = 0; i < pre_requests.size(); i ++) {
		SimpleRequest r1 = pre_requests[i];
		add_edge(2 * i, 2 * i + 1, 1, r1.trip_time);
		for (int j = 0; j < pre_requests.size(); j ++) {
			SimpleRequest r2 = pre_requests[j];
			if (r1.start_time + r1.trip_time + (r1.destination == r2.source ? 0 : dist[r1.destination][r2.source]) <= r2.start_time) {
				add_edge(2 * i + 1, 2 * j, 1, 0);
			}
		}
		for (int j = 0; j < driver_num; j ++) {
			if (drivers[j].time + (r1.source == drivers[j].u ? 0 : dist[drivers[j].u][r1.source]) <= r1.start_time) {
				add_edge(pre_requests.size() * 2 + 2 * j, 2 * i, 1, 0);
			}
			if (r1.start_time + r1.trip_time + (r1.destination == drivers[j].u ? 0 : dist[r1.destination][drivers[j].u]) < TOTAL_TIME_INTERVAL) {
				add_edge(2 * i + 1, pre_requests.size() * 2 + 2 * j + 1, 1, 0);
			}
		}
	}
	for (int i = 0; i < driver_num; i ++) {
		add_edge(source, pre_requests.size() * 2 + i * 2, 1, 0);
		add_edge(pre_requests.size() * 2 + i * 2 + 1, target, 1, 0);
	}
	return max_flow();
}

double offline_value;

void sample_in_distribution() {
	memset(isChosenDriver, 0, sizeof(isChosenDriver));
	driver_num = REAL_DRIVER_NUM;
	for (int i = 0; i < driver_num; i ++) {
		int driver = random_engine() % driver_num;
		while (isChosenDriver[driver]) {
			driver = random_engine() % driver_num;
		}
		isChosenDriver[driver] = true;
		drivers[i].number = driver;
        drivers[i].time = driver_time[driver];
		drivers[i].u = driver_src[driver];
		drivers[i].v = -1;
		drivers[i].length = 0;
		drivers[i].type = 0;
		drivers[i].missions.clear();
		drivers[i].finished_work = 0;
	}
    for (int i = 0; i < regions; i ++)
        for (int j = 0; j < regions; j ++)
            for (int k = 0; k < TOTAL_TIME_INTERVAL; k ++)
				for (int l = 0; l <= driver_num; l ++)
					Prob[i][j][k][l] = initProb[i][j][k][l];
	ond_num = 0;
	pre_num = 0;
	ond_requests.clear();
	pre_requests.clear();
    for (int i = 0; i < regions; i ++)
        for (int j = 0; j < regions; j ++)
            for (int k = 0; k < TOTAL_TIME_INTERVAL; k ++) {
				SimpleRequest sr;
				sr.source = i;
				sr.destination = j;
				sr.start_time = k;
				sr.trip_time = int(V[i][j][k]);
				int cnt = 0;
                double tot_prob = 0;
                for (cnt = driver_num; cnt >= 0; cnt --) {
    				double random_value = 1.0 * (random_engine() % 10000) / 10000;
    				if (random_value < Prob[i][j][k][cnt] - tot_prob) {
                        break;
                    }
                    tot_prob += Prob[i][j][k][cnt];
                }
                if (cnt <= 0) continue;
                while (cnt --) {
					if (pre_ratio <= 1.) {
						double random_value = 1.0 * (random_engine() % 10000) / 10000;
						if (random_value < pre_ratio - eps) {
			                sr.isPre = 1;
		                    sr.order = pre_num ++;
		                    pre_requests.push_back(sr);
			            }
					} else {
		                sr.isPre = 1;
						for (int times = 0; times < int(pre_ratio); times ++){
	                    	sr.order = pre_num ++;
							pre_requests.push_back(sr);
						}
					}
	                sr.isPre = 0;
                    sr.order = ond_num ++;
                    ond_requests.push_back(sr);
				}
			}
    random_shuffle(ond_requests.begin(), ond_requests.end());
    random_shuffle(pre_requests.begin(), pre_requests.end());
    satisfied_ond = satisfied_pre = 0;
	total_pre_value = 0;
	total_ond_value = 0;
	total_satisfied_pre_value = 0;
	total_satisfied_ond_value = 0;
	for (auto tmp : pre_requests) {
		total_pre_value += tmp.trip_time;
	}
	for (auto tmp : ond_requests) {
		total_ond_value += tmp.trip_time;
	}
    if (isOutputMode) {
        printf("We have %d pre-scheduled requests, whose total value is %d.\n", pre_requests.size(), total_pre_value);
        printf("We have %d on-demand requests, whose total value is %d.\n", ond_requests.size(), total_ond_value);
    }
	offline_value = calc_offline();
}

void reset_for_same_distribution() {
	for (int i = 0; i < driver_num; i ++) {
		int driver = drivers[i].number;
        drivers[i].time = driver_time[driver];
		drivers[i].u = driver_src[driver];
		drivers[i].v = -1;
		drivers[i].length = 0;
		drivers[i].type = 0;
		drivers[i].missions.clear();
		drivers[i].finished_work = 0;
	}
    satisfied_ond = satisfied_pre = 0;
	total_satisfied_pre_value = 0;
	total_satisfied_ond_value = 0;
}

void greedy_algorithm() {
	if (isOutputMode) puts("We are using First-Fit algorithm in stage-1.");
    for (int i = 0; i < driver_num; i ++) {
        drivers[i].missions.clear();
        SimpleRequest sr;
        sr.trip_time = 0;
        sr.start_time = drivers[i].time;
        sr.source = sr.destination = drivers[i].u;
        sr.isPre = 1;
        sr.order = -1;
        drivers[i].missions.push_back(sr);
        SimpleRequest sr2;
        sr2.trip_time = 1;
        sr2.start_time = TOTAL_TIME_INTERVAL - 1;
        sr2.source = sr.destination = drivers[i].u;
        sr2.isPre = 1;
        sr2.order = -1;
        drivers[i].missions.push_back(sr2);
    }
    satisfied_pre = 0;
    for (int i = 0; i < pre_num; i ++) {
        for (int j = 0; j < driver_num; j ++) {
			if (drivers[j].time >= pre_requests[i].start_time) continue;
			int pos;
			for (pos = 0; pos < drivers[j].missions.size(); pos ++) {
				if (drivers[j].missions[pos].start_time > pre_requests[i].start_time){
					break;
				}
			}
			if (pos == 0) {
				printRequest(pre_requests[i]);
			}
			assert(pos > 0);
            int arrival_time = drivers[j].missions[pos - 1].start_time + drivers[j].missions[pos - 1].trip_time
				+ calc_time(drivers[j].missions[pos - 1].destination, pre_requests[i].source);
            int ending_time = max(arrival_time, pre_requests[i].start_time) + pre_requests[i].trip_time;
            if (arrival_time <= pre_requests[i].start_time
                && ending_time + calc_time(pre_requests[i].destination, drivers[j].missions[pos].source) <= drivers[j].missions[pos].start_time) {
                drivers[j].missions.insert(drivers[j].missions.begin() + pos, pre_requests[i]);
                satisfied_pre ++;
				total_satisfied_pre_value += pre_requests[i].trip_time;
                break;
            }
        }
    }
	for (int i = 0; i < driver_num; i ++) {
		drivers[i].missions.erase(drivers[i].missions.begin());
	}
	if (isOutputMode) {
    	puts("Performance:--------------------");
        printf("We have handled %d pre-requests of all the %d pre-requests. "
               "Percentage: %.2f%%.\n", satisfied_pre, pre_num, double(100.0 * satisfied_pre / pre_num));
        printf("We have handled %d values of all the %d values in pre-requests. "
               "Percentage: %.2f%%.\n", total_satisfied_pre_value, total_pre_value, double(100.0 * total_satisfied_pre_value / total_pre_value));
   }
}

void heuristic_algorithm() {
    if (isOutputMode) puts("We are using the heuristic algorithm in stage-1.");
    for (int i = 0; i < driver_num; i ++) {
        drivers[i].missions.clear();
        SimpleRequest sr;
        sr.trip_time = 0;
        sr.start_time = drivers[i].time;
        sr.source = sr.destination = drivers[i].u;
        sr.isPre = 1;
        sr.order = -1;
        drivers[i].missions.push_back(sr);
        SimpleRequest sr2;
        sr2.trip_time = 1;
        sr2.start_time = TOTAL_TIME_INTERVAL - 1;
        sr2.source = sr.destination = drivers[i].u;
        sr2.isPre = 1;
        sr2.order = -1;
        drivers[i].missions.push_back(sr2);
    }
    satisfied_pre = 0;
    double bb = 0;
    int cc = 0;
    for (int i = 0; i < pre_num; i ++) {
        if (i % 10 == 0) {
            // printf("We have dealt with %d pre-scheduled requests.\n", i);
        }
		int best_j = -1, best_pos = -1;
		double best_value = 0, best_origin = -(1 << 30);
        for (int j = 0; j < driver_num; j ++) {
			if (drivers[j].time >= pre_requests[i].start_time) continue;
			int pos;
			for (pos = 0; pos < drivers[j].missions.size(); pos ++) {
				if (drivers[j].missions[pos].start_time > pre_requests[i].start_time){
					break;
				}
			}
			assert(pos > 0);
            int arrival_time = drivers[j].missions[pos - 1].start_time + drivers[j].missions[pos - 1].trip_time
				+ calc_time(drivers[j].missions[pos - 1].destination, pre_requests[i].source);
            int ending_time = max(arrival_time, pre_requests[i].start_time) + pre_requests[i].trip_time;
            if (arrival_time <= pre_requests[i].start_time
                && ending_time + calc_time(pre_requests[i].destination, drivers[j].missions[pos].source) <= drivers[j].missions[pos].start_time) {
				int v1 = drivers[j].missions[pos - 1].destination;
				int t1 = drivers[j].missions[pos - 1].start_time + drivers[j].missions[pos - 1].trip_time;
				int v2 = pre_requests[i].destination;
				int t2 = pre_requests[i].start_time + pre_requests[i].trip_time;
				time_stamp ++;
				double tmp1 = pre_computed_table[v1][t1][drivers[j].missions[pos].source][drivers[j].missions[pos].start_time];
                // cout << "---" << v1 << ' ' << t1 << ' ' << drivers[j].missions[pos].source << ' ' << drivers[j].missions[pos].start_time << ' ' << tmp1 << endl;
				double tmp2 = pre_computed_table[v2][t2][drivers[j].missions[pos].source][drivers[j].missions[pos].start_time];
                // cout << v2 << ' ' << t2 << ' ' << drivers[j].missions[pos].source << ' ' << drivers[j].missions[pos].start_time << ' ' << tmp2 << endl;
				tmp2 += pre_requests[i].trip_time;
                // cout << tmp2 << endl;
				time_stamp ++;
				tmp2 += pre_computed_table[v1][t1][pre_requests[i].source][pre_requests[i].start_time];
                // cout << v1 << ' ' << t1 << ' ' << pre_requests[i].source << ' ' << pre_requests[i].start_time << ' ' << tmp2 << endl;
				double value;
				if (Valuemode == "random") {
					value = calc_heuristic_score_function(tmp2 - tmp1);
				} else if (Valuemode == "no") {
					value = tmp2 - tmp1;
				} else {
					assert(false);
				}
				if (best_j == -1 || value > best_value) {
					best_j = j;
					best_pos = pos;
					best_value = value;
                    best_origin = max(best_origin, tmp2 - tmp1);
				}
                // cout << best_value << ' ' << best_origin << endl;
            }
        }
        // best_value += 3.27986;
        if (best_j == -1) continue;
		bool taken = false;
        if (Rejectionmode == "random" && (best_value > 0 || U(random_engine2) < exp(best_value))) taken = true;
        if (Rejectionmode == "strong" && best_value > 0) taken = true;
		if (Rejectionmode == "no") taken = true;
		if (taken) {
			drivers[best_j].missions.insert(drivers[best_j].missions.begin() + best_pos, pre_requests[i]);
			satisfied_pre ++;
			total_satisfied_pre_value += pre_requests[i].trip_time;
		}
        // } // else {
            // cout << best_value << ' ' << best_origin << endl;
        // }
        bb += best_value;
        cc += 1;
    }
    if (isOutputMode) {
        cout << bb / cc << endl;
        cout << VALUE1 / COUNT << endl;
    }
	// necessary operation because of implementation
	for (int i = 0; i < driver_num; i ++) {
		drivers[i].missions.erase(drivers[i].missions.begin());
        for (int j = 0; j < drivers[i].missions.size() - 1; j ++)
            assert(drivers[i].missions[j].start_time <= drivers[i].missions[j + 1].start_time);
	}

	if (isOutputMode) {
    	puts("Performance:--------------------");
        printf("We have handled %d pre-requests of all the %d pre-requests. "
               "Percentage: %.2f%%.\n", satisfied_pre, pre_num, double(100.0 * satisfied_pre / pre_num));
        printf("We have handled %d values of all the %d values in pre-requests. "
               "Percentage: %.2f%%.\n", total_satisfied_pre_value, total_pre_value, double(100.0 * total_satisfied_pre_value / total_pre_value));
   }
}

void GetPlan(const int &v, const int &t, const int &c, const double &x, const SimpleRequest &r) {
    if (v == r.source && t == r.start_time) return;
    vector<pair<int, pair<double, double>>> available_regions;
    for (int d = 0; d < regions; d ++) {
        if ((d == r.source ? 0 : dist[d][r.source]) <= r.start_time - t - dist[v][d]) {
            available_regions.push_back(make_pair(d, make_pair(V[v][d][t], calc_score_value(d, t + dist[v][d], r.source, r.start_time))));
        }
    }
    /*
	cout << t << ' ' << v << ' ' << r.start_time << ' ' << r.source << endl;
	cout << dist[v][r.source] << endl;
	cout << available_regions.size() << endl;
	*/
    assert(v != -1);
    sort(available_regions.begin(), available_regions.end(), cmp_score);
    double d_star = -1;
    int d_index = -1;
    for (auto tmp : available_regions) {
        if (d_index == -1 || d_star < tmp.second.second) {
            d_star = tmp.second.second;
            d_index = tmp.first;
        }
    }
    if (v == r.source && t + 1 <= r.start_time || v != r.source && dist[v][r.source] <= r.start_time - t - 1) {
        double tmp = calc_score_value(v, t + 1, r.source, r.start_time);
        if (d_star < tmp) {
            d_star = tmp;
            d_index = -2;
        }
    }
    int j = -1;
    for (int i = 0; i < available_regions.size(); i ++) {
        if (available_regions[i].second.first + available_regions[i].second.second > d_star) {
            j = i;
        }
    }
    double p = 1.0;
    for (int i = 0; i <= j; i ++) {
        int ai = available_regions[i].first;
        // g[v][ai][t] += x * p * Prob[v][ai][t][1];
        g[v][ai][t] += x * p;
		assert(g[v][ai][t] > -eps && g[v][ai][t] < 1 + eps);
        GetPlan(ai, t + dist[v][ai], c, x * p * Prob[v][ai][t][1], r);
        p *= 1.0 - Prob[v][ai][t][1];
    }
    assert(d_index != -1);
    if (d_index == -2) {
        GetPlan(v, t + 1, c, x * p, r);
    } else {
        GetPlan(d_index, t + dist[v][d_index], c, x * p, r);
    }
}

double TopoValue[REGION_NUM][TOTAL_TIME_INTERVAL * 2];

void TopoUpdate(const int &v, const int &t, const int &c, const double &x, const SimpleRequest &r) {
    if (v == r.source && t == r.start_time) return;
    vector<pair<int, pair<double, double>>> available_regions;
    for (int d = 0; d < regions; d ++) {
        if ((d == r.source ? 0 : dist[d][r.source]) <= r.start_time - t - dist[v][d]) {
            available_regions.push_back(make_pair(d, make_pair(V[v][d][t], calc_score_value(d, t + dist[v][d], r.source, r.start_time))));
        }
    }
    /*
	cout << t << ' ' << v << ' ' << r.start_time << ' ' << r.source << endl;
	cout << dist[v][r.source] << endl;
	cout << available_regions.size() << endl;
	*/
    assert(v != -1);
    sort(available_regions.begin(), available_regions.end(), cmp_score);
    double d_star = -1;
    int d_index = -1;
    for (auto tmp : available_regions) {
        if (d_index == -1 || d_star < tmp.second.second) {
            d_star = tmp.second.second;
            d_index = tmp.first;
        }
    }
    if (v == r.source && t + 1 <= r.start_time || v != r.source && dist[v][r.source] <= r.start_time - t - 1) {
        double tmp = calc_score_value(v, t + 1, r.source, r.start_time);
        if (d_star < tmp) {
            d_star = tmp;
            d_index = -2;
        }
    }
    int j = -1;
    for (int i = 0; i < available_regions.size(); i ++) {
        if (available_regions[i].second.first + available_regions[i].second.second > d_star) {
            j = i;
        }
    }
    double p = 1.0;
    for (int i = 0; i <= j; i ++) {
        int ai = available_regions[i].first;
        // g[v][ai][t] += x * p * Prob[v][ai][t][1];
        g[v][ai][t] += x * p;
        TopoValue[ai][t + dist[v][ai]] += x * p * Prob[v][ai][t][1];
        p *= 1.0 - Prob[v][ai][t][1];
    }
    assert(d_index != -1);
    if (d_index == -2) {
        TopoValue[v][t + 1] += x * p;
    } else {
        TopoValue[d_index][t + dist[v][d_index]] += x * p;
    }
}

priority_queue<Event> que;
vector<vector<SimpleRequest>> ond_requests_time;

void TopoGetPlan(const int &v, const int &t, const int &c, const double &x, const SimpleRequest &r) {
    memset(TopoValue, 0, sizeof(TopoValue));
    TopoValue[v][t] = x;
    for (int tt = t; tt < TOTAL_TIME_INTERVAL; tt ++)
        for (int vv = 0; vv < regions; vv++)
            if (TopoValue[vv][tt] > eps)
                TopoUpdate(vv, tt, c, TopoValue[vv][tt], r);
}

void UpdateDistribution(const int &c, const pair<int, int> &action, const int &t) {
    for (auto sr : offline_requests) {
        g[sr.source][sr.destination][sr.start_time] = 0;
    }
    int va, ta;
    if (action.first) {
        SimpleRequest sr = ond_requests_time[t][action.second];
        va = sr.destination;
        ta = t + sr.trip_time;
        g[sr.source][sr.destination][sr.start_time] = 1;
    } else {
        va = action.second;
        ta = t + (drivers[c].u == va ? 1 : dist[drivers[c].u][va]);
    }
//    cout << va << ' ' << ta << endl;
    SimpleRequest r = drivers[c].missions[drivers[c].finished_work];
//    cout << r.start_time << endl;
//    GetPlan(va, ta, c, 1, r);
    TopoGetPlan(va, ta, c, 1, r);
//    cout << va << ' ' << ta << endl;
    for (auto r : offline_requests) {
        double p = g[r.source][r.destination][r.start_time];
        for (int i = 1; i <= driver_num; i ++) {
            Prob[r.source][r.destination][r.start_time][i] *= (1 - p);
            Prob[r.source][r.destination][r.start_time][i] += p * Prob[r.source][r.destination][r.start_time][i + 1];
			if (Prob[r.source][r.destination][r.start_time][i] < -eps) {
				cout << p << ' ' << Prob[r.source][r.destination][r.start_time][i + 1] << endl;
				assert(false);
			}
        }
    }
}

pair<int, int> ChooseAction(const int &c, const int &v, const int &t) {
    vector<pair<int, int>> actionSet;
    SimpleRequest r = drivers[c].missions[drivers[c].finished_work];
    int best_action = -1;
    double best_value = -1;
	// if (times > 0 && c == 2) {
		// cout << '-' << endl;
	// }
	// if (times > 0 && c == 2) {
	// 	// cout << '-' << endl;
	// }
    for (int i = 0; i < ond_requests_time[t].size(); i ++) {
        if (assigned_driver[i] != -1) continue;
        SimpleRequest sr = ond_requests_time[t][i];
        if (sr.start_time != t || sr.source != v) continue;
        int d = sr.destination;
        if (dist[d][r.source] <= r.start_time - t - sr.trip_time) {
            double tmp = sr.trip_time + calc_score_value(sr.destination, t + sr.trip_time, r.source, r.start_time);
            if (best_action == -1 || tmp > best_value) {
                best_action = i;
                best_value = tmp;
            }
        }
    }
	// if (times > 0 && c == 2) {
	// 	cout << '-' << endl;
	// }
    bool isIlde = false;
	// if (best_action == -1) {
	for (int i = 0; i < regions; i ++) {
        int d = i;
        if ((d == r.source ? 0 : dist[d][r.source]) <= r.start_time - t - (i == v ? 1 : dist[v][i])) {
        // if ((d == r.source ? 0 : dist[d][r.source]) <= r.start_time - t - (i == v ? 1 : dist[v][i])) {
            double tmp = calc_score_value(d, t + (i == v ? 1 : dist[v][i]), r.source, r.start_time);
            if (best_action == -1 || tmp > best_value) {
                best_action = i;
                best_value = tmp;
                isIlde = true;
            }
        }
    }
	// }

// 	if (times > 0 && c == 2) {
// 		cout << '-' << endl;
// 	}
//
//     if (times > 0 && c == 2) {
// 		cout << t << ' ' << v << ' ' << r.start_time << ' ' << r.source << endl;
// 		cout << dist[v][r.source] << endl;
// 		cout << isIlde << endl;
//         cout << r.source << ' ' << dist[drivers[c].u][r.source] << endl;
//         cout << int(isIlde) << ' ' << best_action << endl;
// //		exit(0);
// 	}

    if (isIlde) {
    	drivers[c].u = v;
    	drivers[c].v = best_action;
    	drivers[c].length = (best_action == v ? 1 : dist[v][best_action]) - 1;
        if (drivers[c].length == 0) {
        	drivers[c].u = drivers[c].v;
        	drivers[c].v = -1;
        }
        return make_pair(0, best_action);
    } else {
		// if (c == 0 && DEBUG) {
		// 	SimpleRequest sr = ond_requests_time[t][best_action];
		// 	cerr << sr.trip_time << ' ' << sr.start_time << ' ' << sr.source << ' ' << sr.destination << endl;
		// }
    	satisfied_ond ++;
    	total_satisfied_ond_value += ond_requests_time[t][best_action].trip_time;
    	drivers[c].u = v;
    	drivers[c].v = ond_requests_time[t][best_action].destination;
    	drivers[c].length = ond_requests_time[t][best_action].trip_time - 1;
    	assigned_driver[best_action] = c;
        if (drivers[c].length == 0) {
        	drivers[c].u = drivers[c].v;
        	drivers[c].v = -1;
        }
        return make_pair(1, best_action);
    }
}

void SequentialAlgorithm(const int &t) {
	for (int i = 0; i < regions; i ++)
		for (int j = 0; j < regions; j ++)
			for (int k = 0; k < TOTAL_TIME_INTERVAL; k ++)
				for (int l = 0; l <= driver_num; l ++)
					Prob[i][j][k][l] = initProb[i][j][k][l];
	for (int i = 0; i < ond_requests_time[t].size(); i ++)
		assigned_driver[i] = -1;
	// if (times > 0) cout << t << endl;
    for (int c = 0; c < driver_num; c ++) {
        if (t < drivers[c].time) continue;
        // if (times > 0) cout << c << endl;
        // if (times > 0) cout << '-' << endl;
        time_stamp ++;
        // running on road
        if (drivers[c].length > 0) {
            drivers[c].length --;
            if (drivers[c].length == 0) {
            	drivers[c].u = drivers[c].v;
            	drivers[c].v = -1;
            }
            continue;
        }
        // must pick up a pre-scheduled request
        SimpleRequest r = drivers[c].missions[drivers[c].finished_work];
        if (drivers[c].u == r.source && t == r.start_time) {
            drivers[c].v = r.destination;
            drivers[c].length = r.trip_time - 1;
            drivers[c].finished_work ++;
            if (drivers[c].length == 0) {
            	drivers[c].u = drivers[c].v;
            	drivers[c].v = -1;
            }
			if (c == 0 && DEBUG) {
				// SimpleRequest sr = r;
				// cerr << sr.trip_time << ' ' << sr.start_time << ' ' << sr.source << ' ' << sr.destination << endl;
			}
            continue;
        }
        // if (times > 0) {
		// 	cout << '-' << endl;
		// 	cout << drivers[c].u << ' ' << drivers[c].v << ' ' << drivers[c].length << ' ' << drivers[c].finished_work << ' ' << drivers[c].missions.size() << endl;
		// }
	    // can do ilde drives / can pick up an on-demand request
        pair<int, int> action = ChooseAction(c, drivers[c].u, t);
//        cout << action.first << ' ' << action.second << endl;
        // if (times > 0) cout << '-' << endl;
        UpdateDistribution(c, action, t);
    }
}

void simulation_MDP() {
    if (isOutputMode) puts("We are using MDP algorithm in stage-2.\n");
	// printMissions(0);
	Map.clear();
	ond_requests_time.clear();
    for (int T = 0; T < TOTAL_TIME_INTERVAL; T ++) {
		vector<SimpleRequest> v;
        v.clear();
		ond_requests_time.push_back(v);
	}
	for (auto r : ond_requests) {
		ond_requests_time[r.start_time].push_back(r);
	}
    for (int T = 0; T < TOTAL_TIME_INTERVAL; T ++) {
        if (T > 0 && T % (3600 / TIME_STEP) == 0) printf("Hour : %d\n", T / (3600 / TIME_STEP));
        SequentialAlgorithm(T);
    }
    if (isOutputMode) {
    	puts("Performance:--------------------");
        printf("We have handled %d on-demand requests of all the %d on-demand requests. "
               "Percentage: %.2f%%.\n", satisfied_ond, ond_num, double(100.0 * satisfied_ond / ond_num));
        printf("We have handled %d values of all the %d values in on-demand requests. "
               "Percentage: %.2f%%.\n", total_satisfied_ond_value, total_ond_value, double(100.0 * total_satisfied_ond_value / total_ond_value));
   }
}

int nx, ny;

struct kuhn_munkres {
	int n, w[DRIVER_NUM][DRIVER_NUM], lx[DRIVER_NUM], ly[DRIVER_NUM], m[DRIVER_NUM], way[DRIVER_NUM], sl[DRIVER_NUM];
	bool u[DRIVER_NUM];
	void reset(int x) {
		n = x;
		for (int i = 1; i <= n; i ++)
			for (int j = 1; j <= n; j ++)
				w[i][j] = 0;
	}
	void hungary(int x) {
		m[0] = x;
		int j0 = 0;
		std::fill (sl, sl + n + 1, inf);
		std::fill (u, u + n + 1, false);
		do {
			u[j0] = true;
			int i0 = m[j0], d = inf, j1 = 0;
			for (int j = 1; j <= n; ++ j)
				if (u[j] == false) {
					int cur = -w[i0][j] - lx[i0] - ly[j];
					if (cur < sl[j]) {
						sl[j] = cur;
						way[j] = j0;
					}
					if (sl[j] < d) {
						d = sl[j];
						j1 = j;
					}
				}
			for (int j = 0; j <= n; ++j) {
				if (u[j]) {
					lx[m[j]] += d;
					ly[j] -= d;
				}
				else sl[j] -= d;
			}
			j0 = j1;
		} while (m[j0] != 0);
		do {
			int j1 = way[j0];
			m[j0] = m[j1];
			j0 = j1;
		} while (j0);
	}
	void solve() {
		for (int i = 1; i <= n; ++i) m[i] = lx[i] = ly[i] = way[i] = 0;
		for (int i = 1; i <= n; ++i) hungary (i);
		// int sum = 0;
		// for (int i = 1; i <= n; ++i) sum += w[m[i]][i];
		// return sum;
	}
} KM;

void KM_Algorithm(const int &t) {
	nx = driver_num;
	ny = ond_requests_time[t].size();
	// std::cerr << nx << ' ' << ny << '\n';
	// assert(ond_requests_time[t].size() <= driver_num);
	KM.reset(max(nx, ny));
    for (int c = 0; c < driver_num; c ++) {
        if (t < drivers[c].time) continue;
        // running on road
        if (drivers[c].length > 0) {
            drivers[c].length --;
            if (drivers[c].length == 0) {
            	drivers[c].u = drivers[c].v;
            	drivers[c].v = -1;
            }
			// if (c == 0) {
			// 	cerr << drivers[c].u << ' ' << drivers[c].v << ' ' << drivers[c].length << endl;
			// 	cerr << sr.trip_time << ' ' << sr.start_time << ' ' << sr.source << ' ' << sr.destination << endl;
			// }
            continue;
        }
        // must pick up a pre-scheduled request
        SimpleRequest r = drivers[c].missions[drivers[c].finished_work];
        if (drivers[c].u == r.source && t == r.start_time) {
            drivers[c].v = r.destination;
            drivers[c].length = r.trip_time - 1;
            drivers[c].finished_work ++;
			// if (c == 0) {
			// 	cerr << "driver : " << drivers[c].u << ' ' << drivers[c].v << ' ' << drivers[c].length << endl;
			// 	cerr << "request : " << r.source << ' ' << r.destination << ' ' << r.start_time << ' ' << r.trip_time << "----" << endl;
			// }
            continue;
        }
        // can do ilde drives / can pick up an on-demand request
		bool can_pick = false;
		time_stamp ++;
		for (int i = 0; i < ond_requests_time[t].size(); i ++) {
			SimpleRequest sr = ond_requests_time[t][i];
			if (sr.start_time != t || sr.source != drivers[c].u) continue;
			int d = sr.destination;
			if (dist[d][r.source] <= r.start_time - t - sr.trip_time) {
				if (KMmode == "with future") {
					KM.w[c + 1][i + 1] = sr.trip_time + calc_score_value(sr.destination, t + sr.trip_time, r.source, r.start_time);
				} else if (KMmode == "without future") {
					KM.w[c + 1][i + 1] = sr.trip_time;
				} else {
					assert(false);
				}
				can_pick = true;
			}
		}
		if (can_pick) continue;
		// must go to pick up next pre-scheduled request
        if (drivers[c].u != r.source && t + dist[drivers[c].u][r.source] == r.start_time) {
            drivers[c].v = r.source;
            drivers[c].length = dist[drivers[c].u][r.source] - 1;
	        if (drivers[c].length == 0) {
	        	drivers[c].u = drivers[c].v;
	        	drivers[c].v = -1;
	        }
        }
		// else {
		// 	int d = random_engine() % regions;
		// 	if (d == drivers[c].u) continue;
		// 	while (t + dist[drivers[c].u][d] + (d == r.source ? 0 : dist[d][r.source]) > r.start_time) {
		// 		d = random_engine() % regions;
		// 		if (d == drivers[c].u) break;
		// 	}
		// 	if (d == drivers[c].u) continue;
		// 	drivers[c].v = d;
		// 	drivers[c].length = dist[drivers[c].u][d] - 1;
		// }
    }
	KM.solve();
	// std::cerr << nx << ' ' << ny << '\n';
    for (int i = 1; i <= ny; i ++) {
		int c = KM.m[i];
        if (c > 0 && KM.w[c][i] > 0) {
            c --;
			SimpleRequest sr = ond_requests_time[t][i - 1];
	    	satisfied_ond ++;
	    	total_satisfied_ond_value += sr.trip_time;
	    	drivers[c].u = sr.source;
	    	drivers[c].v = sr.destination;
	    	drivers[c].length = sr.trip_time - 1;
	        if (drivers[c].length == 0) {
	        	drivers[c].u = drivers[c].v;
	        	drivers[c].v = -1;
	        }
			// if (c == 0) {
			// 	cerr << "driver : " << drivers[c].u << ' ' << drivers[c].v << ' ' << drivers[c].length << endl;
			// 	cerr << "request : " << sr.source << ' ' << sr.destination << ' ' << sr.start_time << ' ' << sr.trip_time << ' ' << endl;
			// }
		}
	}
}

void simulation_KM() {
    if (isOutputMode) puts("We are using KM algorithm in stage-2.");
	// for(int i = 0; i < regions; i ++) {
	// 	for (int j = 0; j < regions; j ++)
	// 		cout << dist[i][j] << ' ' ;
	// 	cout << endl;
	// }
	ond_requests_time.clear();
    for (int T = 0; T < TOTAL_TIME_INTERVAL; T ++) {
		vector<SimpleRequest> v;
		ond_requests_time.push_back(v);
	}
	for (auto r : ond_requests) {
		ond_requests_time[r.start_time].push_back(r);
	}
	int max_size = 0;
	for (int T = 0; T < TOTAL_TIME_INTERVAL; T ++) {
		max_size = max(max_size, int(ond_requests_time[T].size()));
	}
    if (isOutputMode) cout << "There are at most " << max_size << " requests in one time step." << endl;
    for (int T = 0; T < TOTAL_TIME_INTERVAL; T ++) {
        //if (T > 0 && T % (3600 / TIME_STEP) == 0) printf("Hour : %d\n", T / (3600 / TIME_STEP));
        KM_Algorithm(T);
    }
    if (isOutputMode) {
    	puts("Performance:--------------------");
        printf("We have handled %d on-demand requests of all the %d on-demand requests. "
               "Percentage: %.2f%%.\n", satisfied_ond, ond_num, double(100.0 * satisfied_ond / ond_num));
        printf("We have handled %d values of all the %d values in on-demand requests. "
               "Percentage: %.2f%%.\n", total_satisfied_ond_value, total_ond_value, double(100.0 * total_satisfied_ond_value / total_ond_value));
   }
}

void arrange_drivers_pre() {
    // any pre-algorithm here
    if (stage1mode == "first-fit") {
		greedy_algorithm();
	} else if (stage1mode == "heuristic") {
	    heuristic_algorithm();
	}
}

void arrange() {
	sample_in_distribution();
    arrange_drivers_pre();
    if (stage2mode == "MDP") {
	    simulation_MDP();
    } else if (stage2mode == "KM") {
	    simulation_KM();
    }
}

double calc_mu() {
	int mx = 0;
	int mn = 1 << 28;
	for (auto tmp : pre_requests) {
		mx = max(mx, tmp.trip_time);
		mn = min(mn, tmp.trip_time);
	}
	return 1.0 * mx / mn;
}

int main() {
    random_engine.seed(seed);
    random_engine2.seed(seed);
    load_data(mode);
	double performance = 0.;
    int avg_pre_num = 0;
    int avg_ond_num = 0;
    int avg_pre_value = 0;
    int avg_ond_value = 0;
	double fbh = -1;
	for (times = 0; times < TOTAL_TIMES; times ++) {
        if (isOutputMode) printf("\nRound %d:\n", times);
		sample_in_distribution();
		performance = 0;
		for (int times2 = 0; times2 < TIMES_FOR_ONE_DISTRIBUTION; times2 ++) {
			reset_for_same_distribution();
		    arrange_drivers_pre();
			double CR = 1.0 * offline_value / total_satisfied_pre_value;
			performance += CR;
	        // avg_pre_num += pre_num;
	        // avg_pre_value += total_pre_value;
			double mu = calc_mu();
			cout << "offline value = " << offline_value << "\tonline value = " << total_satisfied_pre_value << endl;
			cout << "CR = " << CR << " 4 * mu - 1 = " << 4.0 * mu - 1 << endl;
		}
		cout << "Avg CR for this distribution : " << performance / TIMES_FOR_ONE_DISTRIBUTION << endl;
	}
    return 0;
}
// 4 43 34

