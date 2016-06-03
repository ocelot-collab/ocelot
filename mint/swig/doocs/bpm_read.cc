#include <stdio.h>
#include <eq_client.h>

int main(int argc, char* argv[]) {

    char   buf[STRING_LENGTH]; 
    int    rc;
    EqAdr  ea; 
    EqData src;
    EqData dst; 
    EqCall eq;
    
    //char addr[] = "TTF2.DIAG/BPM/1UBC3/CH_X.TD";
    //char addr[] = "TTF2.DIAG/BPM/1UBC3/X";
    char* addr; addr = argv[1];
    // Set the address of a DOOCS channel 
    ea.adr(addr); 

    printf(addr);
    // Make the call
    rc = eq.get(&ea, &src, &dst); 

    // Check the result 
    if (rc) {
        printf("\nRead error %d\n", dst.error());
	}
    else {   
        printf("\nData type %s (%d)\n", dst.type_string(), dst.type());
	printf("Length %d %d\n", dst.length(), dst.array_length());
        printf("\nData is %s\n", dst.get_string(buf, sizeof(buf)));        
    }
    
    double* d_arr = new double[2048];
    //d_arr = dst.get_double_array();
    
    
    for(int i=0; i < dst.length(); i++) {
    	//printf("%f", d_arr[i]);
	d_arr[i] = dst.get_double(i);
	printf("%f\n", d_arr[i]);	
    }
    
    
    //d_arr[0];
    //double x = dst.get_double(0);
    printf("blahblah %f\n", d_arr[0]);	
    delete d_arr;
    
	
    return 0;
}

