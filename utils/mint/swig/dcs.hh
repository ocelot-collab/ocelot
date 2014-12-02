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


class BPM {
public:
	char* id;
	double x;
	double y;
	double s_pos;
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


double get_device_val(char* device_name);

int test_func_1d(Func_1d* f, int n);

Func_1d test_func_1d_2(int n);

MParameters get_parameters();

int get_bpm_val(BPM*);

Orbit get_orbit();

//int track(Field* f, Parameters p);
