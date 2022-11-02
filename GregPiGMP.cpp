
#include <math.h>
#include <stdio.h>
#include <string.h>
//#include <iostream>   //prevents the library attachment
//#include <fstream>
#include <ctime>
//#include <cstring>
#include <gmp.h>


//#define ncells 2000		//remember you need extra cells for accuracy in division using Divlx

/*Developed by Olatunde Adeyemo �2013
GregPiGMP uses Gregory series to calculate Pi as in
 calculate Pi =16 arctan 0.2-4 arctan(1/239)
 * 
 * It uses the GMP arbritrary precision library

let us assume we are working to n dp through int Ndp
*/
//#define USEC_TO_SEC 1.0e-6

//using namespace std;



#include "./fns.cpp"

int main(int argc, char* argv[])
{
	mpf_set_default_prec ( 10000*log(10.0)/log(2.0));
 	mpf_t f,f2,r5,sin,cos,Pi4;
 	mpf_init(f);
 	mpf_init(f2);
 	mpf_init(r5);
 	mpf_init(sin);
 	mpf_init(cos);
 	mpf_init(Pi4);
	char istring[100000];
	
	
	//int Ndp=3000; //  Ndp max ~ncells*6*0.6
	
	mpf_t one, tan, h, t,dup;	
	mpf_init(t);
	mpf_init(one);
	mpf_init(tan);
	mpf_init(h);
	mpf_init(dup);
	
	itoxR(1,&one);

	
	double tsec;

	time_t t1,t2;
        char filename[25]; 	
    
     strcpy(filename, "./PiGregGMP.txt");


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
	fprintf(OutFile,"\n\nPI =\n%s",istring);
	//printf("\n\nPI =\n%s",istring);
	
	
	mpf_neg(dup,r5);
	Divlxr( dup,f,&Pi4);//Cos  Pi4
	//r5.sgn=1;
	//Pi4=r5;	
	
	
	mpf_set(dup,Pi4);	
	CosSin(dup,0,&cos);
	
	crstrXlr(cos,1000,istring);
	fprintf(OutFile,"\ncos(PI/4) =\n");		
	fprintf(OutFile,"%s\n",istring);		
	
	mpf_set(dup,Pi4);
	
	CosSin(dup,1,&sin);		
	//dia =4 ;
	strcpy(istring,"");
	crstrXlr(sin,1000,istring);	
	fprintf(OutFile,"\n\nsin(PI/4) =\n%s\n",istring);	
		
	//sin.sgn *=-1;
	mpf_neg(sin,sin);
	
	addXlr(sin,cos,&f);
	//mpf_out_str (stdout,10,20,f);
	crstrXlr(f,20,istring);
	fprintf(OutFile,"\n\ncos(PI/4)-sin(PI/4) =\n%s",istring);	
    dia =0;  
		
	Divlxr(cos,sin,&tan);
	//tan.sgn=-1;
	mpf_neg(dup,tan);		// sin -ve tan=-1
	addXlr(one,tan,&f);
	//tan.sgn=1;
	addXlr(one,dup,&f2);
	Divlxr(f2,f,&h);
	printf("\n");
	mpf_out_str (stdout,10,20,f);
	printf("\n\n");
	mpf_set(dup,h);
	invtan(dup,&h);
	mpf_set(dup,h);
	addXlr(Pi4,dup,&h);
	itoxR(4,&t);
	mpf_set(dup,h);
	prodxlr(t,dup,&h); 

	time(&t2);
	tsec=difftime(t2, t1);
	
	crstrXlr(h,10000,istring);
	fprintf(OutFile,"\nnew estimate for PI = %s\n",istring);	
      
	fprintf(OutFile,"\n\n total compute time  %lf sec\n",tsec);

}
