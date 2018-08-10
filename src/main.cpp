#include <bits/stdc++.h>

using namespace std;

const int seed = 2333;
const int REGION_NUM = 30;
const int DRIVER_NUM = 40775;

// manual settings
const int DEBUG = 1; // 0 means no
const string mode = "synthetic"; // "file" or "synthetic"
const int MIN_DIST = 100;
const int MAX_DIST = 1000;
const int TOTAL_TIME_INTERVAL = 24 * 60 * 60;
const int MIN_TRIP_TIME = 5 * 60;
const int MAX_TRIP_TIME = 40 * 60;
const int MIN_WAITING_TIME = 5 * 60;
const int MAX_WAITING_TIME = 30 * 60;
const int MIN_ADVANCED_TIME = 20 * 60;
const int MAX_ADVANCED_TIME = 60 * 60;
int driver_num = 40000;
int regions = 21;
int pre_num = 100000;
int ond_num = 200000;

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
};

int dist[REGION_NUM][REGION_NUM];
int driver_time[DRIVER_NUM], driver_src[DRIVER_NUM];
int satisfied_ond = 0, satisfied_pre = 0;

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

vector<SimpleRequest> pre_requests, ond_requests;

int calc_time(const int &l1, const int &l2) {
	return dist[l1][l2];
}

void load_data(const string &mode) {
	if (mode == "file") {

		// TODO

	} else if (mode == "synthetic") {
		for (int i = 0; i < regions; i ++)
			for (int j = 0; j < regions; j ++) {
				if (i <= j) {
					dist[i][j] = rand() % (MAX_DIST - MIN_DIST) + MIN_DIST;
				} else {
					dist[i][j] = dist[j][i];
				}
			}
		for (int i = 0; i < driver_num; i ++) {
			driver_time[i] = rand() % TOTAL_TIME_INTERVAL;
			driver_src[i] = rand() % regions;
			drivers[i] = ActiveDriver(i);
		}
		for (int i = 0; i < pre_num; i ++) {
			SimpleRequest sr;
			sr.order = i;
			sr.start_time = rand() % TOTAL_TIME_INTERVAL; // 00:00 - 23:59
			sr.trip_time = rand() % (MAX_TRIP_TIME - MIN_WAITING_TIME) + MIN_TRIP_TIME; // 5min - 45min
			sr.waiting_time = rand() % (MAX_WAITING_TIME - MIN_WAITING_TIME)+ MIN_WAITING_TIME; // 5min - 30min
			sr.source = rand() % regions;
			sr.destination = rand() % regions;
			sr.isPre = 1;
			pre_requests.push_back(sr);
		}
		for (int i = 0; i < ond_num; i ++) {
			SimpleRequest sr;
			sr.order = i;
			sr.start_time = rand() % TOTAL_TIME_INTERVAL; // 00:00 - 23:59
			sr.trip_time = rand() % (MAX_TRIP_TIME - MIN_WAITING_TIME) + MIN_TRIP_TIME; // 5min - 45min
			sr.waiting_time = rand() % (MAX_WAITING_TIME - MIN_WAITING_TIME)+ MIN_WAITING_TIME; // 5min - 30min
			sr.advanced_time = rand() % (MAX_ADVANCED_TIME - MIN_ADVANCED_TIME)+ MIN_ADVANCED_TIME; // 5min - 10min
			sr.source = rand() % regions;
			sr.destination = rand() % regions;
			sr.isPre = 0;
			ond_requests.push_back(sr);
		}
	}
}

void greedy_algorithm() {
	for (int i = 0; i < driver_num; i ++) {
		drivers[i].missions.clear();
	}
	satisfied_pre = 0;
	for (int i = 0; i < pre_num; i ++) {
		for (int j = 0; j < driver_num; j ++) {
			if (drivers[j].time + calc_time(drivers[j].u, pre_requests[i].source) <= pre_requests[i].start_time + pre_requests[i].waiting_time) { // float operation involved
				drivers[j].missions.push_back(pre_requests[i]);
				drivers[j].time = max(drivers[j].time + calc_time(drivers[j].u, pre_requests[i].source), pre_requests[i].start_time) + pre_requests[i].trip_time;
				drivers[j].u = pre_requests[i].destination;
				satisfied_pre ++;
				break;
			}
		}
	}
	printf("We have handled %d pre-requests of all the %d pre-requests. Percentage: %.2f%%.\n", satisfied_pre, pre_num, double(100.0 * satisfied_pre / pre_num));

	for (int i = 0; i < driver_num; i ++) {
		drivers[i].time = driver_time[i];
		drivers[i].u = driver_src[i];
	}
}

void arrange_drivers_pre() {
	// any pre-algorithm here
	greedy_algorithm();
	// DP();
}

inline bool remainedMission(int driver) {
	return int(drivers[driver].missions.size()) - drivers[driver].finished_work;
}

void printMissions(int driver) {
	printf("---Driver %d has %d missions. Finished %d missions.\n", driver, int(drivers[driver].missions.size()), drivers[driver].finished_work);
	printf("Driver from %d to %d, length: %d\n", drivers[driver].u, drivers[driver].v, drivers[driver].length);
	for (int j = 0; j < drivers[driver].missions.size(); j ++) {
		SimpleRequest sr = drivers[driver].missions[j];
		printf("Mission %d\tRequest %02d\tsource: %d\tdestination: %d\tstart time: %d\toccur time: %d\twaiting time: %d\ttrip time: %d\n", j + 1, sr.order, sr.source, sr.destination, sr.start_time, max(0, sr.start_time - sr.advanced_time), sr.waiting_time, sr.trip_time);
		if (j > 0) printf("dist:%d\n", dist[drivers[driver].missions[j - 1].destination][sr.source]);
	}
	puts("");
}

void printRequest(const SimpleRequest &sr) {
	printf("Request %02d\tsource: %d\tdestination: %d\tstart time: %d\toccur time: %d\twaiting time: %d\ttrip time: %d\n", sr.order, sr.source, sr.destination, sr.start_time, max(0, sr.start_time - sr.advanced_time), sr.waiting_time, sr.trip_time);
	puts("");
}

bool isLegalSequence(const SimpleRequest &request, const vector<SimpleRequest> &tmp2, int position, int driver, const int &T) {
	bool ok;
	if (drivers[driver].type == 0) { // empty car
		int t = max(T, driver_time[driver]);
		int location = (drivers[driver].v != -1 ? drivers[driver].v : drivers[driver].u);
		// if (DEBUG) {
		// 	printf("time:%d\tlocation:%d\n", t, location);
		// }
		if (location != tmp2[0].source) {
			// forward
			ok = true;
			t += dist[location][tmp2[0].source] + (drivers[driver].v != -1 ? drivers[driver].length : 0);
			location = tmp2[0].source;
			for (int i = 0; i < tmp2.size(); i ++) {
				if (t <= tmp2[i].start_time + tmp2[i].waiting_time) {
					t = max(t, tmp2[i].start_time) + tmp2[i].trip_time;
					if (i < tmp2.size() - 1) {
						t += dist[tmp2[i].destination][tmp2[i + 1].source];
					}
				} else {
					ok = false;
					break;
				}
			}
			if (ok) {
				drivers[driver].missions.insert(drivers[driver].missions.begin() + position, request);
				return true;
			}

			if (drivers[driver].v != -1) {
				// backward
				ok = true;
				t = max(T, driver_time[driver]);
				t += dist[drivers[driver].u][drivers[driver].v] - drivers[driver].length;
				if (drivers[driver].u != tmp2[0].source) {
					t += dist[drivers[driver].u][tmp2[0].source];
				}
				location = tmp2[0].source;
				for (int i = 0; i < tmp2.size(); i ++) {
					if (t <= tmp2[i].start_time + tmp2[i].waiting_time) {
						t = max(t, tmp2[i].start_time) + tmp2[i].trip_time;
						if (i < tmp2.size() - 1) {
							t += dist[tmp2[i].destination][tmp2[i + 1].source];
						}
					} else {
						ok = false;
						break;
					}
				}
				if (ok) {
					drivers[driver].missions.insert(drivers[driver].missions.begin() + position, request);
					drivers[driver].length = dist[drivers[driver].u][drivers[driver].v] - drivers[driver].length;
					swap(drivers[driver].u, drivers[driver].v);
					return true;
				}
			}
		} else {
			// already at the start point
			ok = true;
			location = tmp2[0].source;
			t = max(T, driver_time[driver]);
			t += (drivers[driver].v == -1 ? 0 : drivers[driver].length);
			for (int i = 0; i < tmp2.size(); i ++) {
				if (t <= tmp2[i].start_time + tmp2[i].waiting_time) {
					t = max(t, tmp2[i].start_time) + tmp2[i].trip_time;
					if (i < tmp2.size() - 1) {
						t += dist[tmp2[i].destination][tmp2[i + 1].source];
					}
				} else {
					ok = false;
					break;
				}
			}
			if (ok) {
				drivers[driver].missions.insert(drivers[driver].missions.begin() + position, request);
				return true;
			}
		}
	} else { // full car
		int t = max(T, driver_time[driver]);
		int location = tmp2[0].source;
		ok = true;
		t += drivers[driver].length + dist[drivers[driver].v][tmp2[0].source];
		for (int i = 0; i < tmp2.size(); i ++) {
			if (t <= tmp2[i].start_time + tmp2[i].waiting_time) {
				t = max(t, tmp2[i].start_time) + tmp2[i].trip_time;
				if (i < tmp2.size() - 1) {
					t += dist[tmp2[i].destination][tmp2[i + 1].source];
				}
			} else {
				ok = false;
				break;
			}
		}
		if (ok) {
			drivers[driver].missions.insert(drivers[driver].missions.begin() + position, request);
			return true;
		}
	}
	return false;
}

bool assign_to_driver(const SimpleRequest &request, const int &driver, const int &T) {
	vector<SimpleRequest> tmp1, tmp2;
	for (int i = drivers[driver].finished_work; i < drivers[driver].missions.size(); i ++) {
		tmp1.push_back(drivers[driver].missions[i]);
	}
	for (int i = 0; i <= tmp1.size(); i ++) {
		tmp2 = tmp1;
		tmp2.insert(tmp2.begin() + i, request);
		if (isLegalSequence(request, tmp2, drivers[driver].finished_work + i, driver, T)) {
			// if (DEBUG && driver == 23 && request.order == 14171) {
			// 	printf("------------TIME IS %d\n", T);
			// 	printMissions(driver);
			// 	for (int j = 0; j < tmp1.size(); j ++) {
			// 		SimpleRequest sr = tmp1[j];
			// 		printf("Request %02d\tsource: %d\tdestination: %d\tstart time: %d\toccur time: %d\twaiting time: %d\ttrip time: %d\n", sr.order, sr.source, sr.destination, sr.start_time, max(0, sr.start_time - sr.advanced_time), sr.waiting_time, sr.trip_time);
			// 	}
			// }
			return true;
		}
	}
	return false;
}

void arrange_drivers_ond(const int &T, const SimpleRequest &request) {
	for (int i = 0; i < driver_num; i ++) {
		if (assign_to_driver(request, i, T)) return;
	}
}

priority_queue<Event> que;

void simulation() {
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
		e.time = max(0, ond_requests[i].start_time - ond_requests[i].advanced_time);
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
	for (int T = 0; T < 25 * 60 * 60; T ++) {
		if (T > 0 && T % (60 * 60) == 0) printf("Hour : %d\n", T / 60 / 60);
		while (que.top().time == T) {
			Event e = que.top();
			que.pop();
			if (e.type == 0) {
				// handle the orders
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
							if (!(e2.time - 1 <= request.start_time + request.waiting_time)) {
							}
							assert(e2.time - 1 <= request.start_time + request.waiting_time);
							drivers[driver].type = 1;
							drivers[driver].u = request.source;
							drivers[driver].v = request.destination;
							drivers[driver].length = request.trip_time - 1;
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
						printf("\n-----IMPOSSIBLE-----\n");
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
						printf("\n-----IMPOSSIBLE-----\n");
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
						// if (DEBUG && driver == 3015) {
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
	printf("We have handled %d on-demand requests of all the %d on-demand requests. Percentage: %.2f%%.\n", satisfied_ond, ond_num, double(100.0 * satisfied_ond / ond_num));
}

void arrange() {
	arrange_drivers_pre();
	simulation();
}

int main() {
	srand(seed);
	load_data(mode);
	arrange();
	return 0;
}
