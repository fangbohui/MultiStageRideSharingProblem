#include <bits/stdc++.h>

using namespace std;

/*
	In this code:
	driver number: 0-based
	order number: 0-based
*/

int seed = 233;
const int DRIVER_NUM = 40771 + 5;
const int REQUEST_NUM = 214609 + 5;
const int CENTER_NUM = 21;
const long double eps = 1e-9;
long double avg_speed;

long double ABS(long double x) {
	return x > 0 ? x : -x;
}

struct Location {
	long double x, y;
	Location(long double x = 0, long double y = 0) : x(x), y(y) {}
	Location(const Location &location) {
		x = location.x;
		y = location.y;
	}

	long double dist_manh() {
		return ABS(x) + ABS(y);
	}

	long double dist_euc() {
		return sqrt(x * x + y * y);
	}
};

inline Location operator-(const Location &a, const Location &b) {
	return Location(a.x - b.x, a.y - b.y);
}

inline Location operator+(const Location &a, const Location &b) {
	return Location(a.x + b.x, a.y + b.y);
}

inline Location operator/(const Location &a, const long double &b) {
	return Location(a.x / b, a.y / b);
}

inline bool operator==(const Location &a, const Location &b) {
	return ABS(a.x - b.x) < eps && ABS(a.y - b.y) < eps;
}

inline bool operator<(const Location &a, const Location &b) {
	return a.x < b.x - eps || ABS(a.x - b.x) < eps && a.y < b.y - eps;
}

int calc_time(const Location &l1, const Location &l2) {
	return int(1.0 * (l2 - l1).dist_manh() / avg_speed);
}

struct RouteStamp {
	int time;
	Location location;

	void read() {
		double x, y;
		scanf("%d%lf%lf", &time, &x, &y);
		location.x = x;
		location.y = y;
	}
};

struct Request {
	int order, start_time, end_time;
	Location source, destination;

	void read(int number) {
		order = number - 1;
		double sx, sy, dx, dy;
		scanf("%d%d%lf%lf%lf%lf", &start_time, &end_time, &sx, &sy, &dx, &dy);
		source.x = sx;
		source.y = sy;
		destination.x = dx;
		destination.y = dy;
	}
} requests[REQUEST_NUM];

struct Route {
	int driver, order, trip_length;
	vector<RouteStamp> trip_list;

	void read() {
		scanf("%d%d%d", &driver, &order, &trip_length);
		driver --;
		order --;
		for (int i = 0; i < trip_length; i ++) {
			RouteStamp rs;
			rs.read();
			trip_list.push_back(rs);
		}
	}
} routes[REQUEST_NUM];

char requests_filename[23] = "../data/clean_requests";
char routes_filename[21] = "../data/clean_routes";
int route_cnt, request_cnt, drivers, cnt;

void load_data() {
	// read all the routes
	freopen(routes_filename, "r", stdin);
	route_cnt = 0;
	int number = 0;
	while (scanf("%d", &number) == 1) {
		routes[route_cnt ++].read();
		if ((number % 100000) == 0) {
			printf("Read %d routes now.\n", number);
		}
	}
	freopen("/dev/console", "r", stdin);
	printf("Read all the routes. There are %d routes.\n", number);

	// read all the requests
	freopen(requests_filename, "r", stdin);
	request_cnt = 0;
	number = 0;
	while (scanf("%d", &number) == 1) {
		requests[request_cnt ++].read(number);
		if ((number % 100000) == 0) {
			printf("Read %d requests now.\n", number);
		}
	}
	freopen("/dev/console", "r", stdin);
	printf("Read all the routes. There are %d requests.\n", number);
	assert(request_cnt == route_cnt);
	cnt = request_cnt;

	// check the number of drivers
	drivers = -1;
	for (int i = 0; i < cnt; i ++) {
		drivers = max(drivers, routes[i].driver);
	}
	printf("Total drivers: %d\n", drivers);
}

Location driver_src[DRIVER_NUM];
int driver_time[DRIVER_NUM];
vector<Location> points, key_points;

bool check_same(const Location &loc) {
	for (int i = 0; i < key_points.size(); i ++) {
		if (key_points[i] == loc) return true;
	}
	return false;
}

void kmeans(const vector<Location> &points, const int &k) {
	printf("Generating the k-means centers.\n");
	key_points.clear();
	for (int i = 0; i < k; i ++) {
		int ran = rand() % points.size();
		Location loc = points[ran];
		while (check_same(loc)) {
			ran = rand() % points.size();
			loc = points[ran];
		}
		key_points.push_back(loc);
	}
	sort(key_points.begin(), key_points.end());
	vector<Location> loc_zeros;
	vector<int> int_zeros;
	loc_zeros.clear();
	int_zeros.clear();
	for (int i = 0; i < k; i ++) loc_zeros.push_back(Location(0, 0));
	for (int i = 0; i < k; i ++) int_zeros.push_back(0);
	int iter = 0;
	// printf("initial points :\n");
	// for (int i = 0; i < k; i ++)
	// 	printf("%.5f %.5f\n", double(key_points[i].x), double(key_points[i].y));
	vector<Location> means = loc_zeros;
	vector<int> means_cnt = int_zeros;
	while (true) {
		if (iter ++ > 1000) break;
		vector<Location> tmp = key_points;
		means = loc_zeros;
		means_cnt = int_zeros;
		for (int i = 0; i < points.size(); i ++) {
			int min_point = -1;
			for (int j = 0; j < k; j ++) {
				if (min_point == -1 || (points[i] - key_points[min_point]).dist_euc() > (points[i] - key_points[j]).dist_euc()) {
					min_point = j;
				}
			}
			means[min_point] = means[min_point] + points[i];
			means_cnt[min_point] += 1;
		}
		for (int i = 0; i < k; i ++) {
			key_points[i] = means[i] / means_cnt[i];
		}
		// for (int i = 0; i < k; i ++) printf("%.5f %.5f\n", double(key_points[i].x), double(key_points[i].y));

		// check the change
		sort(key_points.begin(), key_points.end());
		bool isSame = true;
		for (int i = 0; i < k; i ++) {
			if (!(ABS(key_points[i].x - tmp[i].x) < eps && ABS(key_points[i].y - tmp[i].y) < eps)) {
				isSame = false;
				break;
			}
		}
		if (isSame) break;
	}
	for (int i = 0; i < k; i ++) {
		for (int j = i + 1; j < k; j ++) {
			if (means_cnt[j] > means_cnt[i]) {
				swap(means_cnt[i], means_cnt[j]);
				swap(key_points[i], key_points[j]);
			}
		}
	}
	printf("Total k-means iteration : %d\n", iter);
	for (int i = 0; i < k; i ++) {
		printf("Key point %d:\t%.6f\t%.6f\tclustered points: %d\n", i, double(key_points[i].x), double(key_points[i].y), means_cnt[i]);
	}
}

void do_statistics() {
	// check the average distance and time
	freopen("distance.txt", "w", stdout);
	long double dist = 0;
	long long time = 0;
	for (int i = 0; i < cnt; i ++) {
		long double tmp1 = 0;
		long double tmp2 = routes[i].trip_list[routes[i].trip_length - 1].time - routes[i].trip_list[0].time;
		for (int j = 1; j < routes[i].trip_length; j ++) {
			tmp1 += (routes[i].trip_list[j].location - routes[i].trip_list[j - 1].location).dist_manh();
		}
		dist += tmp1;
		time += tmp2;
	}
	freopen("/dev/tty", "w", stdout);
	long double avg_dist = 1.0 * dist / cnt;
	long double avg_time = 1.0 * time / cnt;
	avg_speed = avg_dist / avg_time;
	printf("Average distance per route is %.8f.\n", double(avg_dist));
	printf("Average time per route is %.3f seconds.\n", double(avg_time));
	printf("Average speed per second is %.8f.\n", double(avg_speed));

	// check the start points of drivers
	freopen("start_points.txt", "w", stdout);
	for (int i = 0; i < drivers; i ++) {
		driver_time[i] = 48 * 60 * 60;
	}
	for (int i = 0; i < cnt; i ++) {
		if (routes[i].trip_list[0].time < driver_time[routes[i].driver]) {
			driver_time[routes[i].driver] = routes[i].trip_list[0].time;
			driver_src[routes[i].driver] = routes[i].trip_list[0].location;
		}
	}
	freopen("/dev/tty", "w", stdout);
	printf("Calculated all the start points of all the drivers.\n");
	printf("\n");

	// k-means
	points.clear();
	for (int i = 0; i < cnt; i ++) {
		points.push_back(requests[i].source);
		points.push_back(requests[i].destination);
	}
	kmeans(points, CENTER_NUM);
}

int calc_center(const Location &loc) {
	int belonged = -1;
	for (int i = 0; i < CENTER_NUM; i ++) {
		if (belonged == -1 || (loc - key_points[belonged]).dist_euc() > (loc - key_points[i]).dist_euc()) {
			belonged = i;
		}
	}
	return belonged;
}

struct SimpleRequest {
	int order, start_time, end_time, waiting_time, trip_time;
	int source, destination;

	void assign(Request &req, Route &route) {
		order = req.order;
		start_time = req.start_time;
		end_time = req.end_time;
		trip_time = route.trip_list[route.trip_length - 1].time - route.trip_list[0].time;
		waiting_time = end_time - start_time - trip_time;
		assert(waiting_time >= 0);
		source = calc_center(req.source);
		destination = calc_center(req.destination);
	}
};

inline bool operator<(const SimpleRequest &sr1, const SimpleRequest &sr2) {
	return sr1.start_time + sr1.waiting_time < sr2.start_time + sr2.waiting_time;
}

struct ActiveDriver {
	int driver, time;
	Location location;
	ActiveDriver(int driver = -1) : driver(driver), time(driver_time[driver]), location(driver_src[driver]) {}
};

vector<ActiveDriver> valid_drivers;
vector<SimpleRequest> pre_requests, ond_requests; // ond means on-demand ...
int pre_num, ond_num, valid_driver_num;
vector<int> region_cnt[CENTER_NUM][CENTER_NUM];

void sampling(long double driver_prop, long double pre_prop, long double ond_prop) {
	// choose the drivers
	for (int i = 0; i < drivers; i ++) {
		if (1.0 * (rand() % 10000) / 10000 < driver_prop) {
			valid_drivers.push_back(ActiveDriver(i));
		}
	}
	valid_driver_num = valid_drivers.size();
	printf("We have sampled %d drivers.\n", valid_driver_num);

	// choose the pre-requests and on-demand-requests
	pre_requests.clear();
	ond_requests.clear();
	for (int i = 0; i < cnt; i ++) {
		long double random_value = 1.0 * (rand() % 10000) / 10000;
		if (random_value < pre_prop - eps) {
			SimpleRequest sr;
			sr.assign(requests[i], routes[i]);
			pre_requests.push_back(sr);
			region_cnt[sr.source][sr.destination].push_back(sr.trip_time);
		} else if (random_value > pre_prop - eps && random_value < pre_prop + ond_prop - eps) {
			SimpleRequest sr;
			sr.assign(requests[i], routes[i]);
			ond_requests.push_back(sr);
			region_cnt[sr.source][sr.destination].push_back(sr.trip_time);
		}
	}
	pre_num = pre_requests.size();
	ond_num = ond_requests.size();
	printf("We have sampled %d pre-requests.\n", pre_num);
	printf("We have sampled %d on-demand requests.\n", ond_num);

	double aa = 0, bb = 0;
	for (int i = 0; i < 21; i ++) {
		for (int j = 0; j < 21; j ++) {
			int tmp_cnt = 0;
			int tmp_cnt2 = 0;
			for (int k = 0; k < region_cnt[i][j].size(); k ++) {
				if (region_cnt[i][j][k] >= 180) {
					tmp_cnt += region_cnt[i][j][k];
					tmp_cnt2 ++;
				}
			}
			printf("%.5f\t", 1.0 * tmp_cnt / tmp_cnt2);
			aa += tmp_cnt;
			bb += tmp_cnt2;
		}
		printf("\n");
	}
	cout << aa / bb << endl;

	for (int i = 0; i < CENTER_NUM; i ++) {
		for (int j = 0; j < CENTER_NUM; j ++) {
			printf("%d\t", region_cnt[i][j].size());
		}
		printf("\n");
	}

    freopen("../data/driver.txt", "w", stdout);
    printf("%d\n", valid_driver_num);
    for (int i = 0; i < valid_driver_num; i ++) {
        printf("%d %d\n", valid_drivers[i].time, calc_center(valid_drivers[i].location));
    }
	freopen("/dev/tty", "w", stdout);
    printf("Output %d drivers to file.\n", drivers);

    freopen("../data/pre_requests.txt", "w", stdout);
    printf("%d\n", pre_num);
    for (int i = 0; i < pre_num; i ++) {
        printf("%d %d %d %d %d\n", pre_requests[i].start_time, pre_requests[i].end_time, pre_requests[i].waiting_time, pre_requests[i].source, pre_requests[i].destination);
    }
	freopen("/dev/tty", "w", stdout);
    printf("Output %d pre-requests to file.\n", pre_num);

    freopen("../data/ond_requests.txt", "w", stdout);
    printf("%d\n", ond_num);
    for (int i = 0; i < ond_num; i ++) {
        printf("%d %d %d %d %d\n", ond_requests[i].start_time, ond_requests[i].end_time, ond_requests[i].waiting_time, ond_requests[i].source, ond_requests[i].destination);
    }
	freopen("/dev/tty", "w", stdout);
    printf("Output %d on-demand requests to file.\n", ond_num);
}

int main() {
	srand(seed);
	load_data();
	do_statistics();
	sampling(1, 0.2, 0.8); // driver, pre, on-demand
}
