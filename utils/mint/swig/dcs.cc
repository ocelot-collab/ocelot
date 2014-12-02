/* File: example.c */

#include "dcs.hh"

double pi = 3.14;

string version() { return "0.0"; }
char* info() { return "info"; }

MParameters::MParameters(char* desc) { this->energy = 0.0; }

BPM::BPM(char* id) {this->id = id; this->x = 0.0; this->y = 0.0;}
char* BPM::info() { return id;}

Orbit::Orbit() {}
BPM Orbit::get(char* id) {
	BPM b("none");
	//return this->bpms[0];
	for(vector<BPM>::size_type i = 0; i != bpms.size(); i++) {
	    /* std::cout << someVector[i]; ... */
		if ( !strcmp(bpms[i].id, id) ) return bpms[i];
	}
	return b;
}

Func_1d::Func_1d(int n) {this->n = n; f = new double[n]; for(int i=0; i<n;i++) f[i] = 0; }
double Func_1d::get(int i){ return f[i]; }
double* Func_1d::get(){ return f; }
double Func_1d::set(int i, double val){ f[i] = val; }
double Func_1d::sum(){
		double s = 0;
		for(int i=0;i<n;i++) s+= f[i];
		return s;
	}


int test_func_1d(Func_1d* g, int n){
	for (int i=0; i<n;i++) g->f[i] = sin(i/10.0);
	return 0.0;
}

Func_1d test_func_1d_2(int n){
	Func_1d f2(n);
	for (int i=0; i<n;i++) f2.f[i] = cos(i/2.0);
	return f2;
}

MParameters get_parameters() {
	MParameters p("flash");
	p.energy = 678;
	return p;
}

double get_device_val(char* device_name) {
	cout << "debug: getting device value for " << device_name << endl;
	return -12.0;
}

int get_bpm_val(BPM* b) {
	b->x = rand() % 100;
	b->y = rand() % 100;
}

Orbit get_orbit() {
	Orbit orb;

	BPM b1("bpm1"); get_bpm_val(&b1);//b1.x = 1.0;
	BPM b2("bpm2"); get_bpm_val(&b2);//b2.x = 2.0;

	orb.bpms.push_back(b1);
	orb.bpms.push_back(b2);

	return orb;
}

//int track(Field *f, Parameters p) { f->f[0] += p.x; }
