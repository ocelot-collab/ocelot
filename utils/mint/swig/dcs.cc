/* simple c++ interface to doocs for python binding */

#include "dcs.hh"
#include <eq_client.h>

bool debug = false;

double pi = 3.14;

string version() { return "0.0"; }
char* info() { return "info"; }

MParameters::MParameters(char* desc) { this->energy = 0.0; }

BPM::BPM(char* id) {
	this->id = (char*) malloc(strlen(id) + 1); 
	memcpy(this->id, id, strlen(id));
	this->channel_x = (char*) malloc( strlen(id) + strlen("/X.FLASH1") + 1);
	memcpy(this->channel_x, id, strlen(id));
	memcpy(this->channel_x + strlen(id), "/X.FLASH1", strlen("/X.FLASH1") + 1);//+1 to copy the null-terminator
	this->channel_y = (char*) malloc( strlen(id) + strlen("/Y.FLASH1") + 1);
	memcpy(this->channel_y, id, strlen(id));
	memcpy(this->channel_y + strlen(id), "/Y.FLASH1", strlen("/Y.FLASH1") + 1);//+1 to copy the null-terminator
	this->channel_z = (char*) malloc( strlen(id) + strlen("/Z_POS") + 1);
	memcpy(this->channel_z, id, strlen(id));
	memcpy(this->channel_z + strlen(id), "/Z_POS", strlen("/Z_POS") + 1);//+1 to copy the null-terminator

	this->x = 0.0; 
	this->y = 0.0;
	this->z_pos = 0.0;
}

Device::Device(char* id) {
	this->id = (char*) malloc(strlen(id) + 1); 
	memcpy(this->id, id, strlen(id));
	this->channel_z = (char*) malloc( strlen(id) + strlen("/Z_POS") + 1);
	memcpy(this->channel_z, id, strlen(id));
	memcpy(this->channel_z + strlen(id), "/Z_POS", strlen("/Z_POS") + 1);//+1 to copy the null-terminator

	this->z_pos = 0.0;
}


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

Func_1d get_device_td(char* device_name) {

	char buf[STRING_LENGTH]; 
	int rc;
    EqAdr  ea; 
    EqData src;
    EqData dst; 
    EqCall eq;
    
    char* addr; addr = device_name;
    ea.adr(addr); 

    if (debug) printf(addr);
    rc = eq.get(&ea, &src, &dst); 

    if (rc) {
      printf("\nRead error %d\n", dst.error());
	}
    else {   
        if (debug) printf("\nData type %s (%d)\n", dst.type_string(), dst.type());
	if (debug) printf("Length %d %d\n", dst.length(), dst.array_length());
	if (debug)  printf("\nData is %s\n", dst.get_string(buf, sizeof(buf)));        
    }
    
	int n = dst.length();
	Func_1d f2(n);
	for (int i=0; i<n;i++) f2.f[i] = dst.get_double(i);
	return f2;
        
    

}

MParameters get_parameters() {
	MParameters p("flash");
	p.energy = 678;
	return p;
}

double get_device_val(char* device_name) {
  if (debug) cout << "debug: getting device value for " << device_name << endl;


  char buf[STRING_LENGTH]; 
  int rc;
  EqAdr  ea; 
  EqData src;
  EqData dst; 
  EqCall eq;
  
  char* addr; addr = device_name;
  ea.adr(addr); 
  
  rc = eq.get(&ea, &src, &dst); 

  if (rc) {
    printf("\nRead error %d\n", dst.error());
    return 0.0;
  }
  else {  return dst.get_double();     }
    
}

int set_device_val(char* device_name, double val) {

  EqData  src;
  int    rc;
  EqAdr  ea; 
  
  rc = src.set(val);

  if (debug) cout << "debug: set_device_val " << device_name << " : " << rc << endl;

  cout << src.get_double() << endl;
  cout<< src.type_string() << endl;

  EqData dst; 
  EqCall eq;
    
  char* addr; addr = device_name;
  printf(addr);
  ea.adr(addr); 
  //rc = eq.get(&ea, &src, &dst); 
  rc = eq.set(&ea, &src, &dst); 
  
  if (rc) {
    printf("\nWrite error %d\n", dst.error());
    return 1;
  }
  
  return 0;
}

int get_bpm_val(BPM* b) {
  if (debug) cout << "debug: getting bpm reading for " << b->id << endl;
  char   buf[STRING_LENGTH]; 
  int    rc;
  EqAdr  ea; 
  EqData src;
  EqData dst; 
  EqCall eq;
  
  char* addr; 
  
  addr = b->channel_x;
  ea.adr(addr); 
  if (debug) cout << "device address : " << addr << endl;
  rc = eq.get(&ea, &src, &dst); 
  if (rc) { printf("\nRead error %d\n", dst.error()); }
  else {  b->x = dst.get_double(); }
  
  addr = b->channel_y;
  ea.adr(addr); 
  if (debug) cout << "device address : " << addr << endl;
  rc = eq.get(&ea, &src, &dst); 
  if (rc) { printf("\nRead error %d\n", dst.error()); }
  else {  b->y = dst.get_double(); }
  
  addr = b->channel_z;
  ea.adr(addr); 
  if (debug) cout << "device address : " << addr << endl;
  rc = eq.get(&ea, &src, &dst); 
  if (rc) { printf("\nRead error %d\n", dst.error()); }
  else {  b->z_pos = dst.get_double(); }
  
	
}

int get_device_info(Device* d) {
  if (debug) cout << "debug: getting device info for " << d->id << endl;
  char   buf[STRING_LENGTH]; 
  int    rc;
  EqAdr  ea; 
  EqData src;
  EqData dst; 
  EqCall eq;
  
  char* addr; 	
  addr = d->channel_z;
  ea.adr(addr); 
  if (debug) cout << "device address : " << addr << endl;
  rc = eq.get(&ea, &src, &dst); 
  if (rc) { printf("\nRead error %d\n", dst.error()); }
  else {  d->z_pos = dst.get_double(); }

}

Orbit get_orbit() {
	Orbit orb;

	BPM b1("bpm1"); get_bpm_val(&b1);//b1.x = 1.0;
	BPM b2("bpm2"); get_bpm_val(&b2);//b2.x = 2.0;

	orb.bpms.push_back(b1);
	orb.bpms.push_back(b2);

	return orb;
}
