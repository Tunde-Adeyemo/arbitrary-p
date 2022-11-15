
#include <math.h>
#include <stdio.h>
#include <string.h>
//#include <iostream>
#include <fstream>
#include <time.h>
#include <cstring>
#include <gmp.h>



/*Developed by Olatunde Adeyemo �2013
GregPi uses Gregory series to calculate Pi as in
 calculate Pi =16 arctan 0.2-4 arctan(1/239)

let us assume we are working to n dp through int Ndp
*/
//#define USEC_TO_SEC 1.0e-6

using namespace std;



#include "./fns.cpp"

int main(int argc, char* argv[])
{
	int Ndp = 100000;
	mpf_set_default_prec ( Ndp*log(10.0)/log(2.0));
 	mpf_t f,f2,r5,sin,cos,Pi4;
 	mpf_init(f);
 	mpf_init(f2);
 	mpf_init(r5);
 	mpf_init(sin);
 	mpf_init(cos);
 	mpf_init(Pi4);
 	
	char istring[100000];
	
		
	mpf_t one, tan, h, t,dup;	
	mpf_init(one);
 	mpf_init(tan);
 	mpf_init(h);
 	mpf_init(t);
 	mpf_init(dup);
	itoxR(1,&one);

	
	double tsec;

	time_t t1,t2;
        char filename[25]; 	
    
     strcpy(filename, ".//PiGreg.txt");


    FILE* OutFile = fopen(filename, (char *) &"w");
    if (OutFile == NULL) {
        printf("Cannot create output files\n");
        printf("Usage: Execute native/seq_demo2 from ..."); 
        exit(0);
    }
    
	time(&t1);
	
	fprintf(OutFile,"\n\n  calculate Pi =16 arctan 0.2-4 arctan(1/239) \n%s",istring);


	itoxR(5,&dup);
	Divlxr( dup,one,&r5);
	invtan(r5,&f);
	
	crstrXlr(f,200,istring);
	printf("\n\ninvtan 0.2\n%s\n",istring);
	
	itoxR(16,&r5);
	mpf_set(dup,f);
	prodxlr(r5,dup,&f);

	itoxR(239,&dup);
	Divlxr( dup,one,&r5);
	
	crstrXlr(r5,200,istring);
	printf("\n\n 1/239\n%s\n",istring);
	
	invtan(r5,&f2);
	
	crstrXlr(f2,200,istring);
	printf("\n\ninvtan 1/239\n%s\n",istring);
	
	itoxR(-4,&r5);
	//r5.sgn=-1;
	mpf_set(dup,f2);
	prodxlr(r5,dup,&f2);
	mpf_set(dup,f);		
	addXlr(dup,f2,&f);		//	*** Pi
			
	crstrXlr(f,10000,istring);
	fprintf(OutFile,"\r\n\nPI =\n%s",istring);

	mpf_neg(dup,r5);	
	Divlxr( dup,f,&Pi4);//Cos  Pi4	
	//r5.sgn=1;
	//Pi4=r5;	
	mpf_set(dup,Pi4);
	
	CosSin(dup,0,&cos);
	crstrXlr(cos,1000,istring);
	fprintf(OutFile,"\r\n\ncos(PI/4) =\n%s",istring);
	mpf_set(dup,Pi4);
	CosSin(dup,1,&sin);		
	crstrXlr(sin,1000,istring);
	fprintf(OutFile,"\n\nsin(PI/4) =\n%s",istring);	
	
	//sin.sgn *=-1;
	mpf_neg(dup,sin);
	
	addXlr(dup,cos,&f);
	crstrXlr(f,20,istring);
	fprintf(OutFile,"\n\ncos(PI/4)-sin(PI/4) =\n%s",istring);
    	
	
	Divlxr(cos,sin,&tan);
	//tan.sgn=-1;
	mpf_neg(dup,tan);
	addXlr(one,dup,&f);	
	//tan.sgn=1;
	addXlr(one,tan,&f2);
	
	Divlxr(f2,f,&h);
	mpf_set(dup,h);
	
	invtan(dup,&h);
	mpf_set(dup,h);
	addXlr(Pi4,dup,&h);
	itoxR(4,&t);
	mpf_set(dup,h);
	prodxlr(t,h,&h); 

	time(&t2);
	tsec=difftime(t2, t1);	
	
	crstrXlr(h,10000,istring);
	fprintf(OutFile,"\r\nnew estimate for PI =\r%s",istring);
      
	fprintf(OutFile,"\r\n\n total compute time  %lf sec\n",tsec);

}
