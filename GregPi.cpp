
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <time.h>
//#include <cstring>


#define ncells 2000		//remember you need extra cells for accuracy in division using Divlx

/*Developed by Olatunde Adeyemo �2013
GregPi uses Gregory series to calculate Pi as in
 calculate Pi =16 arctan 0.2-4 arctan(1/239)

let us assume we are working to n dp through int Ndp
*/
//#define USEC_TO_SEC 1.0e-6

using namespace std;



struct NxlongR
{
	int mant[ncells]; //each cell of array mant carries 6 digits of our extended real number
	int exp;	//exp is the exponent of our first non zero digit of mant[1]
	int sgn;	//overall sign -1 or 1 or if  NxlongR==0 sgn =0
};

ifstream InFile;   // used with  rdXlr()

#include "./fnct.ext"

int main(int argc, char* argv[])
{
 	NxlongR f,f2,r5,sin,cos,Pi4;
	char istring[8000];
	
	int Ndp=2000; //  Ndp max ~ncells*6*0.6
	
	NxlongR one, tan, h, t;	
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


	itoxR(5,&r5);
	Divlxr( r5,one,&r5,Ndp);
	invtan(r5,&f,Ndp);
	
	crstrXlr(f,200,istring);
	printf("\n\ninvtan 0.2\n%s\n",istring);
	
	itoxR(16,&r5);
	prodxlr(r5,f,&f,Ndp);

	itoxR(239,&r5);
	Divlxr( r5,one,&r5,Ndp);
	
	crstrXlr(r5,200,istring);
	printf("\n\n 1/239\n%s\n",istring);
	
	invtan(r5,&f2,Ndp);
	
	crstrXlr(f2,200,istring);
	printf("\n\ninvtan 1/239\n%s\n",istring);
	
	itoxR(4,&r5);
	r5.sgn=-1;
	prodxlr(r5,f2,&f2,Ndp);		
	addXlr(f,f2,&f,Ndp);		//	*** Pi
			
	crstrXlr(f,Ndp,istring);
	fprintf(OutFile,"\r\n\nPI =\n%s",istring);

	Divlxr( r5,f,&r5,Ndp);//Cos  Pi4
	r5.sgn=1;
	Pi4=r5;	
	
	CosSin(Pi4,0,Ndp,&cos);
	crstrXlr(cos,Ndp,istring);
	fprintf(OutFile,"\r\n\ncos(PI/4) =\n%s",istring);
	
	CosSin(Pi4,1,Ndp,&sin);		
	crstrXlr(sin,Ndp,istring);
	fprintf(OutFile,"\n\nsin(PI/4) =\n%s",istring);
	
	sin.sgn *=-1;
	
	addXlr(sin,cos,&f,Ndp);
	crstrXlr(f,20,istring);
	fprintf(OutFile,"\n\ncos(PI/4)-sin(PI/4) =\n%s",istring);
       
	
	
	Divlxr(cos,sin,&tan,Ndp);
	tan.sgn=-1;
	addXlr(one,tan,&f,Ndp);
	tan.sgn=1;
	addXlr(one,tan,&f2,Ndp);
	Divlxr(f2,f,&h,Ndp);
	invtan(h,&h,Ndp);
	addXlr(Pi4,h,&h,Ndp);
	itoxR(4,&t);
	prodxlr(t,h,&h,Ndp); 

	time(&t2);
	tsec=difftime(t2, t1);
	
	crstrXlr(h,Ndp,istring);
	fprintf(OutFile,"\r\nnew estimate for PI =\r%s",istring);
      
	fprintf(OutFile,"\r\n\n total compute time  %lf sec\n",tsec);

}
