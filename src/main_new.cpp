#include <bits/stdc++.h>
#include <random>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>

using namespace std;

default_random_engine random_engine;
const double eps = 1e-9;
const int seed = 2333;
const int REGION_NUM = 30;
const int DRIVER_NUM = 105;
const int REQUEST_NUM = 400000;
const int TIME_STEP = 60;
const int TOTAL_TIME_INTERVAL = 24 * (3600 / TIME_STEP); // TODO make the timestep be 1s when available

const string mode = "file"; // "file" or "synthetic"
const string stage1mode = "heuristic"; // "first-fit" or "heuristic"
const string stage2mode = "MDP";

// manual settings for "synthetic" mode
const int DEBUG = 1; // 0 means no
const int MIN_DIST = 100;
const int MAX_DIST = 1000;
const int MIN_TRIP_TIME = 5 * 60;
const int MAX_TRIP_TIME = 40 * 60;
//const int MIN_WAITING_TIME = 5 * 60;
//const int MAX_WAITING_TIME = 30 * 60;
int driver_num = 100;
int regions = 21;
int pre_num = 100000;
int ond_num = 200000;

// manual settings for "file" mode
double pre_ratio = 0.05; // 0.05;
double ond_ratio = 0.95; // 0.95;
double driver_ratio = 1.;
int days = 10;

// heuristic settings
const double a = 1.;
const double b = 1.;
const double k = 1.;
uniform_real_distribution<double> U(0.0, 1.0);

double VALUE1 = 0;
int COUNT = 0;

double calc_heuristic_score_function(double delta) {
	double alpha = U(random_engine);
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
int driver_time[DRIVER_NUM], driver_src[DRIVER_NUM];
int satisfied_ond = 0, satisfied_pre = 0;
int Count[REGION_NUM][REGION_NUM][TOTAL_TIME_INTERVAL * 2];
int CountValue[REGION_NUM][REGION_NUM][TOTAL_TIME_INTERVAL * 2];
double Prob[REGION_NUM][REGION_NUM][TOTAL_TIME_INTERVAL * 2][DRIVER_NUM];
double V[REGION_NUM][REGION_NUM][TOTAL_TIME_INTERVAL * 2];

struct ActiveDriver {
    int driver, time;
    int u, v, length, finished_work;
    int type; // 0 means ilde, 1 means full
    vector<SimpleRequest> missions;
    ActiveDriver(int driver = -1) : driver(driver), time(driver_time[driver]), u(driver_src[driver]), v(-1), length(-1), finished_work(0), type(0), missions() {}
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
        scanf("%d", &driver_num);
        driver_num = 100; // TODO
        for (int i = 0; i < driver_num; i ++) {
            scanf("%d%d", &driver_time[i], &driver_src[i]);
            driver_time[i] = 0;
            drivers[i] = ActiveDriver(i);
        }
        freopen("../data/all_the_requests.txt", "r", stdin);
        pre_num = 0;
        ond_num = 0;
        int total_num;
        scanf("%d", &total_num);
        printf("%d\n", total_num);
        for (int i = 0; i < total_num; i ++) {
            SimpleRequest sr;
            scanf("%d%d%d%d%d", &sr.start_time, &sr.trip_time, &sr.waiting_time, &sr.source, &sr.destination);
            if (sr.trip_time <= 180) continue;
            sr.start_time /= TIME_STEP;
            sr.trip_time /= TIME_STEP;
//            sr.trip_time = dist[sr.source][sr.destination];
            sr.waiting_time /= TIME_STEP;
            sr.waiting_time = 0;
            long double random_value = 1.0 * (random_engine() % 10000) / 10000;
            int day = random_engine() % days;
            if (random_value < pre_ratio - eps) {
                sr.isPre = 1;
                ond_requests_days[day].push_back(sr);
                if (day == days - 1) {
                    sr.order = pre_num ++;
                    pre_requests.push_back(sr);
                } else {
                    offline_requests.push_back(sr);
                }
            } else if (random_value < pre_ratio + ond_ratio - eps) {
                sr.isPre = 0;
                ond_requests_days[day].push_back(sr);
                if (day == days - 1) {
                    sr.order = ond_num ++;
                    ond_requests.push_back(sr);
                } else {
                    offline_requests.push_back(sr);
                }
            }
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
            driver_time[i] = 0;
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
    for (int day = 0; day < days - 1; day ++) {
        for (int i = 0; i < ond_requests_days[day].size(); i++) {
            SimpleRequest sr = ond_requests_days[day][i];
            Count[sr.source][sr.destination][sr.start_time] ++;
        }
        for (int i = 0; i < regions; i ++)
            for (int j = 0; j < regions; j ++)
                for (int k = 0; k < TOTAL_TIME_INTERVAL; k ++)
                    if (Count[i][j][k] >= 1) {
                        for (int l = 1; l <= Count[i][j][k]; l++)
                            Prob[i][j][k][l] += 1.0 / days;
                        Count[i][j][k] = 0;
                    }
    }
    for (int day = 0; day < days - 1; day ++) {
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
    printf("We have %d pre-scheduled requests.\n", pre_requests.size());
    printf("We have %d on-demand requests.\n", ond_requests.size());
    printf("We have %d offline requests.\n", offline_requests.size());
}
int total_pre_value = 0, total_ond_value = 0;
int total_satisfied_pre_value = 0, total_satisfied_ond_value = 0;

void greedy_algorithm() {
	puts("We are using First-Fit algorithm in stage-1.");
    for (int i = 0; i < driver_num; i ++) {
        drivers[i].missions.clear();
        SimpleRequest sr;
        sr.trip_time = 0;
        sr.start_time = 0;
        sr.source = sr.destination = driver_src[i];
        sr.isPre = 1;
        sr.order = -1;
        drivers[i].missions.push_back(sr);
        SimpleRequest sr2;
        sr2.trip_time = 1;
        sr2.start_time = TOTAL_TIME_INTERVAL - 1;
        sr2.source = sr.destination = driver_src[i];
        sr2.isPre = 1;
        sr2.order = -1;
        drivers[i].missions.push_back(sr2);
    }
    satisfied_pre = 0;
    for (int i = 0; i < pre_num; i ++) {
        total_pre_value += pre_requests[i].trip_time;
        for (int j = 0; j < driver_num; j ++) {
			int pos;
			for (pos = 0; pos < drivers[j].missions.size(); pos ++) {
				if (drivers[j].missions[pos].start_time >= pre_requests[i].start_time){
					break;
				}
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
    printf("We have handled %d pre-requests of all the %d pre-requests. "
           "Percentage: %.2f%%.\n", satisfied_pre, pre_num, double(100.0 * satisfied_pre / pre_num));
    printf("We have handled %d values of all the %d values in pre-requests. "
           "Percentage: %.2f%%.\n", total_satisfied_pre_value, total_pre_value, double(100.0 * total_satisfied_pre_value / total_pre_value));
}

inline bool cmp_score(const pair<int, pair<double, double>> &a, const pair<int, pair<double, double>> &b) {
    return a.second.first + a.second.second > b.second.first + b.second.second;
}

double g[REGION_NUM][REGION_NUM][TOTAL_TIME_INTERVAL];

map<pair<int, int>, pair<int, double>> Map;
int time_stamp = 0;

double calc_score_value(const int &v, const int &t, const SimpleRequest &request) {
//	if (check) cout << v << ' ' << t << endl;
    if (v == request.source && t == request.start_time) return 0;
    if (Map[make_pair(v, t)].first == time_stamp) {
        return Map[make_pair(v, t)].second;
    }
    vector<pair<int, pair<double, double>>> available_regions;
    for (int d = 0; d < regions; d ++) {
        if ((d == request.source ? 0 : dist[d][request.source]) <= request.start_time - t - dist[v][d]) {
            available_regions.push_back(make_pair(d, make_pair(V[v][d][t], calc_score_value(d, t + dist[v][d], request))));
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
    if (v == request.source && t + 1 <= request.start_time) {
        d_star = max(d_star, calc_score_value(v, t + 1, request));
    } else if (v != request.source && dist[v][request.source] <= request.start_time - t - 1) {
        d_star = max(d_star, calc_score_value(v, t + 1, request));
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
        F += p * Prob[v][ai][t][1] * (V[v][ai][t] + calc_score_value(ai, t + dist[v][ai], request));
        p *= 1.0 - Prob[v][ai][t][1];
    }
    F += p * d_star;
    Map[make_pair(v, t)].first = time_stamp;
    Map[make_pair(v, t)].second = F;
    return F;
}

void heuristic_algorithm() {
	puts("We are using the heuristic algorithm in stage-1.");
    for (int i = 0; i < driver_num; i ++) {
        drivers[i].missions.clear();
        SimpleRequest sr;
        sr.trip_time = 0;
        sr.start_time = 0;
        sr.source = sr.destination = driver_src[i];
        sr.isPre = 1;
        sr.order = -1;
        drivers[i].missions.push_back(sr);
        SimpleRequest sr2;
        sr2.trip_time = 1;
        sr2.start_time = TOTAL_TIME_INTERVAL - 1;
        sr2.source = sr.destination = driver_src[i];
        sr2.isPre = 1;
        sr2.order = -1;
        drivers[i].missions.push_back(sr2);
    }
    satisfied_pre = 0;
    double bb = 0;
    int cc = 0;
    for (int i = 0; i < pre_num; i ++) {
        if (i > 0 && i % 1 == 0) {
            printf("We have dealt with %d pre-scheduled requests.\n", i);
        }
        total_pre_value += pre_requests[i].trip_time;
		int best_j = -1, best_pos = -1;
		double best_value = 0, best_origin = -(1 << 30);
        for (int j = 0; j < driver_num; j ++) {
			int pos;
			for (pos = 0; pos < drivers[j].missions.size(); pos ++) {
				if (drivers[j].missions[pos].start_time >= pre_requests[i].start_time){
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
				double tmp1 = calc_score_value(v1, t1, drivers[j].missions[pos]);
                // cout << "---" << v1 << ' ' << t1 << ' ' << drivers[j].missions[pos].source << ' ' << drivers[j].missions[pos].start_time << ' ' << tmp1 << endl;
				double tmp2 = calc_score_value(v2, t2, drivers[j].missions[pos]);
                // cout << v2 << ' ' << t2 << ' ' << drivers[j].missions[pos].source << ' ' << drivers[j].missions[pos].start_time << ' ' << tmp2 << endl;
				tmp2 += pre_requests[i].trip_time;
                // cout << tmp2 << endl;
				time_stamp ++;
				tmp2 += calc_score_value(v1, t1, pre_requests[i]);
                // cout << v1 << ' ' << t1 << ' ' << pre_requests[i].source << ' ' << pre_requests[i].start_time << ' ' << tmp2 << endl;
				double value = calc_heuristic_score_function(tmp2 - tmp1);
				if (best_j == -1 || value > best_value) {
					best_j = j;
					best_pos = pos;
					best_value = value;
                    best_origin = max(best_origin, tmp2 - tmp1);
				}
                // cout << best_value << ' ' << best_origin << endl;
            }
        }
        best_value += 7.46297;
        if (best_j == -1) continue;
        if (best_value > 0 || U(random_engine) < exp(best_value)){
    		drivers[best_j].missions.insert(drivers[best_j].missions.begin() + best_pos, pre_requests[i]);
    		satisfied_pre ++;
    		total_satisfied_pre_value += pre_requests[i].trip_time;
        } else {
            // cout << best_value << ' ' << best_origin << endl;
        }
        bb += best_value;
        cc += 1;
    }
    // cout << bb / cc << endl;
    // cout << VALUE1 / COUNT << endl;
	// necessary operation because of implementation
	for (int i = 0; i < driver_num; i ++) {
		drivers[i].missions.erase(drivers[i].missions.begin());
        for (int j = 0; j < drivers[i].missions.size() - 1; j ++)
            assert(drivers[i].missions[j].start_time <= drivers[i].missions[j + 1].start_time);
	}
    printf("We have handled %d pre-requests of all the %d pre-requests. "
           "Percentage: %.2f%%.\n", satisfied_pre, pre_num, double(100.0 * satisfied_pre / pre_num));
    printf("We have handled %d values of all the %d values in pre-requests. "
           "Percentage: %.2f%%.\n", total_satisfied_pre_value, total_pre_value, double(100.0 * total_satisfied_pre_value / total_pre_value));
}

void GetPlan(const int &v, const int &t, const int &c, const double &x, const SimpleRequest &r) {
    if (v == r.source && t == r.start_time) return;
    vector<pair<int, pair<double, double>>> available_regions;
    for (int d = 0; d < regions; d ++) {
        if ((d == r.source ? 0 : dist[d][r.source]) <= r.start_time - t - dist[v][d]) {
            available_regions.push_back(make_pair(d, make_pair(V[v][d][t], calc_score_value(d, t + dist[v][d], r))));
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
        double tmp = calc_score_value(v, t + 1, r);
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
        g[v][ai][t] += x * p * Prob[v][ai][t][1];
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
            available_regions.push_back(make_pair(d, make_pair(V[v][d][t], calc_score_value(d, t + dist[v][d], r))));
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
        double tmp = calc_score_value(v, t + 1, r);
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
        g[v][ai][t] += x * p * Prob[v][ai][t][1];
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
        SimpleRequest sr = ond_requests[action.second];
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
        double p = g[r.source][r.destination][r.start_time] / Prob[r.source][r.destination][r.start_time][1];
        for (int i = 1; i <= driver_num; i ++) {
            Prob[r.source][r.destination][r.start_time][i] *= (1 - p);
            Prob[r.source][r.destination][r.start_time][i] += p * Prob[r.source][r.destination][r.start_time][i + 1];
        }
    }
}

pair<int, int> ChooseAction(const int &c, const int &v, const int &t) {
    vector<pair<int, int>> actionSet;
    SimpleRequest r = drivers[c].missions[drivers[c].finished_work];
    int best_action = -1;
    double best_value = -1;
    for (int i = 0; i < ond_requests.size(); i ++) {
        if (assigned_driver[i] != -1) continue;
        SimpleRequest sr = ond_requests[i];
        if (sr.start_time != t || sr.source != v) continue;
        int d = sr.destination;
        if (dist[d][r.source] <= r.start_time - t - sr.trip_time) {
            double tmp = sr.trip_time + calc_score_value(sr.destination, t + sr.trip_time, r);
            if (best_action == -1 || tmp > best_value) {
                best_action = i;
                best_value = tmp;
            }
        }
    }
    bool isIlde = false;
    for (int i = 0; i < regions; i ++) {
        int d = i;
        if ((d == r.source ? 0 : dist[d][r.source]) <= r.start_time - t - (i == v ? 1 : dist[v][i])) {
            double tmp = calc_score_value(d, t + (i == v ? 1 : dist[v][i]), r);
            if (best_action == -1 || tmp > best_value) {
                best_action = i;
                best_value = tmp;
                isIlde = true;
            }
        }
    }
    /*
    if (c == 96) {
		cout << t << ' ' << v << ' ' << r.start_time << ' ' << r.source << endl;
		cout << dist[v][r.source] << endl;
		cout << isIlde << endl;
        cout << r.source << ' ' << dist[drivers[c].u][r.source] << endl;
        cout << int(isIlde) << ' ' << best_action << endl;
//		exit(0);
	}
	*/
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
    	satisfied_ond ++;
    	total_satisfied_ond_value += ond_requests[best_action].trip_time;
    	drivers[c].u = v;
    	drivers[c].v = ond_requests[best_action].destination;
    	drivers[c].length = ond_requests[best_action].trip_time - 1;
    	assigned_driver[best_action] = c;
        if (drivers[c].length == 0) {
        	drivers[c].u = drivers[c].v;
        	drivers[c].v = -1;
        }
        return make_pair(1, best_action);
    }
}

void SequentialAlgorithm(const int &t) {
    for (int c = 0; c < driver_num; c ++) {
//        cout << c << endl;
//        cout << '-' << endl;
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
            continue;
        }
//        cout << '-' << endl;
        // can do ilde drives / can pick up an on-demand request
        pair<int, int> action = ChooseAction(c, drivers[c].u, t);
//        cout << action.first << ' ' << action.second << endl;
//        cout << '-' << endl;
        UpdateDistribution(c, action, t);
    }
}
priority_queue<Event> que;

void simulation_MDP() {
	puts("We are using MDP algorithm in stage-2.\n");
    for (int T = 0; T < TOTAL_TIME_INTERVAL; T ++) {
        if (T > 0 && T % (3600 / TIME_STEP) == 0) printf("Hour : %d\n", T / (3600 / TIME_STEP));
        SequentialAlgorithm(T);
    }
    printf("We have handled %d on-demand requests of all the %d on-demand requests. "
           "Percentage: %.2f%%.\n", satisfied_ond, ond_num, double(100.0 * satisfied_ond / ond_num));
    printf("We have handled %d values of all the %d values in on-demand requests. "
           "Percentage: %.2f%%.\n", total_satisfied_ond_value, total_ond_value, double(100.0 * total_satisfied_ond_value / total_ond_value));
}

bool isLegalSequence_impatient(const SimpleRequest &request, const vector<SimpleRequest> &tmp2, int driver, const int &T) {
    if (drivers[driver].type == 0 && drivers[driver].v == -1 && drivers[driver].u == tmp2[0].source) { // empty car
        // already at the start point
        bool ok = true;
        int t = max(T, driver_time[driver]);
        for (int i = 0; i < tmp2.size(); i ++) {
            if (t <= tmp2[i].start_time) {
                t = max(t, tmp2[i].start_time) + tmp2[i].trip_time;
                if (i < tmp2.size() - 1) {
                    t += dist[tmp2[i].destination][tmp2[i + 1].source];
                }
            } else {
                ok = false;
                assert(i <= 1);
                break;
            }
        }
        if (ok) {
			total_satisfied_ond_value += request.trip_time;
			satisfied_ond ++;
            drivers[driver].missions.insert(drivers[driver].missions.begin() + drivers[driver].finished_work, request);
            return true;
        }
    }
    return false;
}

bool assign_to_driver_impatient(const SimpleRequest &request, const int &driver, const int &T) {
    vector<SimpleRequest> tmp;
    for (int i = drivers[driver].finished_work; i < drivers[driver].missions.size(); i ++) {
        tmp.push_back(drivers[driver].missions[i]);
    }
    tmp.insert(tmp.begin(), request);
    return isLegalSequence_impatient(request, tmp, driver, T);
}

void arrange_drivers_pre() {
    // any pre-algorithm here
    if (stage1mode == "first-fit") {
		greedy_algorithm();
	} else if (stage1mode == "heuristic") {
	    heuristic_algorithm();
	}
}

void arrange_drivers_ond(const int &T, const SimpleRequest &request) {
    for (int i = 0; i < driver_num; i ++) {
        if (assign_to_driver_impatient(request, i, T)) return;
    }
}

void simulation() {
	puts("We are using normal algorithm in stage-2.\n");
    while (!que.empty()) que.pop();
    for (int i = 0; i < driver_num; i ++) {
        Event e;
        e.type = 1;
        e.number = i;
        e.time = driver_time[i];
        que.push(e);
        drivers[i].type = 0;
        drivers[i].u = driver_src[i];
        drivers[i].v = -1;
        drivers[i].length = -1;
    }
    // ond_num = 50;
    for (int i = 0; i < ond_num; i ++) {
        Event e;
        e.type = 0;
        e.number = i;
        e.time = ond_requests[i].start_time;
        que.push(e);
    }

    if (DEBUG) {
        // printf("Distance map:\n");
        // for (int i = 0; i < regions; i ++) {
        // 	for (int j = 0; j < regions; j ++) {
        // 		printf("%d\t", dist[i][j]);
        // 	}
        // 	printf("\n");
        // }

        // for (int i = 0; i < driver_num; i ++) {
        // 	printf("Driver %d: mission number: %d\n", i, int(drivers[i].missions.size()));
        // 	printf("start_region %d: start_time: %d\n", driver_src[i], driver_time[i]);
        // 	for (int j = 0; j < drivers[i].missions.size(); j ++) {
        // 		SimpleRequest sr = drivers[i].missions[j];
        // 		printf("source: %d\tdestination: %d\tstart time: %d\twaiting time: %d\ttrip time: %d\n", sr.source, sr.destination, sr.start_time, sr.waiting_time, sr.trip_time);
        // 	}
        // }

        // for (int i = 0; i < ond_num; i ++) {
        // 	SimpleRequest sr = ond_requests[i];
        // 	printf("On Demand Request %02d: \n", i);
        // 	printf("source: %d\tdestination: %d\tstart time: %d\twaiting time: %d\ttrip time: %d\tadvanced time: %d\n", sr.source, sr.destination, sr.start_time, sr.waiting_time, sr.trip_time, sr.advanced_time);
        // }
    }

    satisfied_ond = 0;
    // in one second:
    //		0.orders occur and assign them
    //		1.drivers move 1 sec OR stay still
    for (int T = 0; T < TOTAL_TIME_INTERVAL; T ++) {
        if (T > 0 && T % (3600 / TIME_STEP) == 0) printf("Hour : %d\n", T / TIME_STEP);
        while (que.top().time == T) {
            Event e = que.top();
            que.pop();
            if (e.type == 0) {
                // handle the orders
                // if (e.number == ) printf("%d\n", e.number);
                arrange_drivers_ond(T, ond_requests[e.number]);
            } else if (e.type == 1) {
                // handle the drivers
                int driver = e.number;
                if (remainedMission(driver) == 0) {
                    // No work for now. Don't move.
                    assert(drivers[driver].v == -1 && drivers[driver].length == -1);
                    Event e2 = e;
                    e2.time ++;
                    que.push(e2);
                } else if (drivers[driver].v != -1 && drivers[driver].length > 1) {
                    // running on road
                    Event e2 = e;
                    e2.time ++;
                    drivers[driver].length --;
                    que.push(e2);
                } else if (drivers[driver].v != -1 && drivers[driver].length == 1) { // A driver is about to finish one's trip towards a source/destination.
                    Event e2 = e;
                    e2.time ++;
                    // pick up the passenger
                    SimpleRequest request;
                    if (remainedMission(driver) != 0) {
                        request = drivers[driver].missions[drivers[driver].finished_work];
                    }
                    if (drivers[driver].type == 0 && drivers[driver].v == request.source) {
                        if (e2.time < request.start_time) {
                            // early arrival
                            drivers[driver].type = 0;
                            drivers[driver].u = request.source;
                            drivers[driver].v = -1;
                            drivers[driver].length = -1;
                        } else {
                            // pick up successfully
                            // if (DEBUG && driver == 3015) {
                            // 	printf("%d\n", driver_time[driver]);
                            // 	printf("Hurry Pick-up driver %d Time %d Request %d\n", driver, T, request.order);
                            // 	printMissions(driver);
                            // 	printf("-------source: %d\tdestination: %d\tstart time: %d\twaiting time: %d\ttrip time: %d\n", request.source, request.destination, request.start_time, request.waiting_time, request.trip_time);
                            // }
                            // if (!(e2.time - 1 <= request.start_time + request.waiting_time)) {
                            // }
                            assert(e2.time - 1 <= request.start_time);
                            drivers[driver].type = 1;
                            drivers[driver].u = request.source;
                            drivers[driver].v = request.destination;
                            drivers[driver].length = request.trip_time - 1;
							if (drivers[driver].length == 0) {
	                            drivers[driver].type = 0;
	                            drivers[driver].u = request.destination;
	                            drivers[driver].v = -1;
								drivers[driver].length = -1;
							}
                        }
                    } else if (drivers[driver].type == 0 && drivers[driver].v != request.source) {
                        // change the target
                        drivers[driver].type = 0;
                        drivers[driver].u = drivers[driver].v;
                        drivers[driver].v = request.source;
                        drivers[driver].length = dist[drivers[driver].u][drivers[driver].v];
                    } else if (drivers[driver].type == 1) {
                        // drop off the passenger
                        drivers[driver].type = 0;
                        drivers[driver].finished_work ++;
                        if (remainedMission(driver) == 0) {
                            // no more work now, stay in v
                            drivers[driver].u = drivers[driver].v;
                            drivers[driver].v = -1;
                            drivers[driver].length = -1;
                        } else {
                            // has other missions
                            // TODO maybe change the logic here
                            if (drivers[driver].v != drivers[driver].missions[drivers[driver].finished_work].source) {
                                drivers[driver].u = drivers[driver].v;
                                drivers[driver].v = drivers[driver].missions[drivers[driver].finished_work].source;
                                drivers[driver].length = calc_time(drivers[driver].u, drivers[driver].v); // ???
                            } else {
                                drivers[driver].u = drivers[driver].v;
                                drivers[driver].v = -1;
                                drivers[driver].length = -1;
                            }
                            // if (driver == 29 && drivers[driver].missions[drivers[driver].finished_work - 1].order == 217) {
                            // 	printf("-----------------I AM HERE. TIME IS %d\n", T);
                            // 	printf("-----------------I AM HERE. u IS %d\n", drivers[driver].u);
                            // 	printf("-----------------I AM HERE. v IS %d\n", drivers[driver].v);
                            // 	printf("-----------------I AM HERE. length IS %d\n", drivers[driver].length);
                            // }
                        }
                        if (!(request.isPre)) satisfied_ond ++;
                    } else {
                        printf("\n-----IMPOSSIBLE1-----\n");
                    }
                    que.push(e2);
                } else if (drivers[driver].v != -1 && drivers[driver].length == 0) {
                    // Special case : only happens in such case that a driver reverses when he just left
                    Event e2 = e;
                    e2.time ++;
                    SimpleRequest request = drivers[driver].missions[drivers[driver].finished_work];
                    if (drivers[driver].type == 0 && drivers[driver].v == request.source) { // pick up the passenger TODO
                        drivers[driver].type = 1;
                        drivers[driver].u = request.source;
                        drivers[driver].v = request.destination;
                        drivers[driver].length = request.trip_time - 1;
                    } else if (drivers[driver].type == 0 && drivers[driver].v != request.source) { // change the target
                        drivers[driver].type = 0;
                        drivers[driver].u = drivers[driver].v;
                        drivers[driver].v = request.source;
                        drivers[driver].length = dist[drivers[driver].u][drivers[driver].v] - 1;
                    } else {
						assert(false);
                        if (DEBUG) {
                            printf("Waiting Pick-up driver %d Time %d Request %d\n", driver, T, request.order);
                            printMissions(driver);
                            printf("-------source: %d\tdestination: %d\tstart time: %d\twaiting time: %d\ttrip time: %d\n", request.source, request.destination, request.start_time, request.waiting_time, request.trip_time);
                        }
                        printf("\n-----IMPOSSIBLE2-----\n");
						assert(false);
                    }
                    que.push(e2);
                } else if (drivers[driver].v == -1 && drivers[driver].length == -1) {
                    // A car is waiting at u
                    Event e2 = e;
                    e2.time ++;
                    assert(drivers[driver].type == 0 && remainedMission(driver) != 0);
                    SimpleRequest request = drivers[driver].missions[drivers[driver].finished_work];
                    if (drivers[driver].u != request.source) {
                        drivers[driver].type = 0;
                        drivers[driver].v = request.source;
                        drivers[driver].length = dist[drivers[driver].u][request.source] - 1;
                    } else if (e2.time - 1 >= drivers[driver].missions[drivers[driver].finished_work].start_time) {
                        // if (DEBUG && driver == 25) {
                        // 	printf("Waiting Pick-up driver %d Time %d Request %d\n", driver, T, request.order);
                        // 	printMissions(driver);
                        // 	printf("-------source: %d\tdestination: %d\tstart time: %d\twaiting time: %d\ttrip time: %d\n", request.source, request.destination, request.start_time, request.waiting_time, request.trip_time);
                        // }
                        assert(e2.time - 1 <= request.start_time + request.waiting_time);
                        drivers[driver].type = 1;
                        drivers[driver].u = request.source;
                        drivers[driver].v = request.destination;
                        drivers[driver].length = request.trip_time - 1;
                    }
                    que.push(e2);
                }
            }
        }
    }
    printf("We have handled %d on-demand requests of all the %d on-demand requests. "
           "Percentage: %.2f%%.\n", satisfied_ond, ond_num, double(100.0 * satisfied_ond / ond_num));
    printf("We have handled %d values of all the %d values in on-demand requests. "
           "Percentage: %.2f%%.\n", total_satisfied_ond_value, total_ond_value, double(100.0 * total_satisfied_ond_value / total_ond_value));
}

void arrange() {
	for (int i = 0; i < ond_num; i ++) {
		total_ond_value += ond_requests[i].trip_time;
	}
    arrange_drivers_pre();
    if (stage2mode == "normal") {
	    simulation(); // only for file now; need to fix later
	} else if (stage2mode == "MDP") {
	    simulation_MDP();
    }
}

int main() {
    random_engine.seed(seed);
    load_data(mode);
    arrange();
    return 0;
}
