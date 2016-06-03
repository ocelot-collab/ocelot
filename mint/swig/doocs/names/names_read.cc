#include <stdio.h>
#include <eq_client.h>

int main(int argc, char *argv[]) {
    char * addr;
    char   buf[STRING_LENGTH]; 
    int    rc;
    EqAdr  ea; 
    EqData src;
    EqData dst; 
    EqCall eq;

    int data_type;
    float f1;
    float f2;
    time_t tm;
    char * com;
    int index, count;

    if (argc > 1)
	{
		addr = argv[1];
	}
	else {
		//addr = "TEST.CAMERA/TTF2ICCD3/IEE1394.CAM1/*";
	        addr = "TEST.CAMERA/*";
	}

    ea.adr(addr);

    rc = eq.names(&ea, &dst);

    if (rc) 
        printf("Error %d\n", dst.error());
    else
        {
            count = dst.length();
            printf("total count is %d\n", count);
            for (index=0; index < count ; index++)
            {
                dst.get_ustr(&data_type, &f1, &f2, &tm, &com, index);
                printf("%d %s %d %d\n", index, com, data_type, sizeof(com));
            }
        }
    return 0;
}

