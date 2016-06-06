/* File: example.h */

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstring>

using namespace std;

string version();

char* info();

extern double pi;

class MParameters {
public:
	MParameters(char * desc);
	double energy;
};


class Device {
public:
	char* id;
	char* channel_z;
	double z_pos;
public:
	Device(char*);
};

class BPM {
public:
	char* id;
	char* channel_x;
	char* channel_y;
	char* channel_z;
	double x;
	double y;
	double z_pos;
public:
	BPM(char*);
	char* info();
};


class Orbit {
public:
	Orbit();
	vector<BPM> bpms;
	BPM get(char*);
};


class Func_1d {
public:
	Func_1d(int n);
	double* f;
	int n;
	double* get();
	double get(int i);
	double set(int i, double v);
	double sum();
};


int get_device_info(Device*);

double get_device_val(char* device_name);
int set_device_val(char* device_name, double val);

Func_1d get_device_td(char* device_name);

int test_func_1d(Func_1d* f, int n);

Func_1d test_func_1d_2(int n);



MParameters get_parameters();

int get_bpm_val(BPM*);

Orbit get_orbit();

//int track(Field* f, Parameters p);
