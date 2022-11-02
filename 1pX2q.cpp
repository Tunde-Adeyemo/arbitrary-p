

#include <math.h>
#include <cstdio>
#include <string.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <cstring>



#define ncells 2000		//remember you need extra cells for accuracy in division using Divlx

/*	
 * 	1pX2q.cpp
 * Developed by Olatunde Adeyemo ©2022
 * 
 * let us assume we are working to n dp through int Ndp
*/
//#define USEC_TO_SEC 1.0e-6

using namespace std;



struct NxlongR
{
	int mant[ncells]; //each cell of array mant carries 6 digits of our extended real number
	int exp;	//exp is the exponent of our first non zero digit of mant[1]
	int sgn;	//overall sign -1 or 1 or if  NxlongR==0 sgn =0
};

//ifstream InFile;   // used with  rdXlr()

#include "./fnct.ext"
/*
 * (1+U)^q = 1+sigmai((u^i . productj(q-j+1)/(i!))
 * 
 * the normal implementation is  u = -+x^2 
 * here we assume  0<x<1 giving
 * 
 *  int((1+x^2)^q )dx= x + sigmai((x^2i+1 . productj(q-j+1)/(i! . (2i+1))
 *  int((1-x^2)^q )dx= x + sigmai((-x^2)i+1 . productj(q-j+1)/(i! . (2i+1))
 * 
 * or with x >1, 
 * (1+U)^q = (U(1+U^-1)) ^q = x^2q *(1+-x^-2)^q
 * int((1+x^2)^q )dx= x^(2q+1)/(2q+1)+sigmai(x^(2(q-i)+1) . productj(2q-j+1)/(i! . 2(q-i)+1)
 * 
 * int((1-x^2)^q )dx= x^(1-2q)/(2q-1) +sigmai((-1^i)x^(2qi+1) . productj(j)/(i! . (2qi+1)))
 * 
 */
 
int main(int argc, char* argv[])
{
	//FILE* OutFile;
 	NxlongR one,X,x2,q,B,fct,prd,k2,termx,tIx,smx,smIx,bufi,Kx2;
	char istring[8000], pmstr[10];
	
	int Ndp=1000, xlt1=0,pm,ep,expo;
	long a,b;
	double Xdbl;
	
	double tsec;

	time_t t1,t2;
        char filename[25]; 	
    
     strcpy(filename, ".//1pX2q.txt");


    FILE* OutFile = fopen(filename, (char *) &"w");
    if (OutFile == NULL) {
        printf("Cannot create output files\n");
        printf("Usage: Execute native/seq_demo2 from ..."); 
        exit(0);
    }
    
	time(&t1);
	
	
	itoxR(1,&one);
    
    /*
     * got to input x, q
     * q read as a/b , a,b integers
     */
     printf("\nenter value x  in the form of 2 entries,\
	\nan upto 15 digit +ve mantissa\n\
	an integer exponent\n\
	enter mantissa eg 2.592\n");
	
    scanf(" %la",&Xdbl);
    printf("\nenter exponent eg -23 or 7\n");
	scanf("%i",&expo);
	//printf("\nXdbl in %lf",Xdbl);
	X.sgn = 1;
	if (Xdbl <0)
	{
		X.sgn =-1;
		Xdbl*= -1;
	}
		
	while(Xdbl<1.0)
	{
		Xdbl*=10;
		ep--;
	}

	while(Xdbl>=10.0)
	{
		Xdbl/=10;
		ep++;
	}
	//printf("\nXdbl %lf",Xdbl);
	expo+=ep;
	Xdbl*=1e5;
	//printf("\nXdbl *1e5 %lf",Xdbl);
	X.mant[1]=Xdbl;
	Xdbl-=X.mant[1];
	Xdbl*=1e6;
	X.mant[2]=Xdbl;
	Xdbl-=X.mant[2];
	Xdbl*=1e6;
	X.mant[3]=Xdbl;
	X.exp=expo;
	
	if (expo >= 0)	
		xlt1=1;
			
	//printf("\nX[1] %d\n",X.mant[1]);
	/*crstrXlr(X,60,istring);
	printf("\n X in\n%s",istring);*/
	
	printf("\nenter rational value q  in the form of 2 integer entries a and b\n \
	 q = a/b  b may be 1,\nenter a\n");
	 scanf("%ld",&a);
	 printf("enter b\n");
	 scanf("%ld",&b);
	 itoxR(a,&q);
	 itoxR(b,&B);
	 Divlxr( B,q,&q,Ndp);        // ** "real" q used in (1-u)^q	
	 
	 
	printf("\nenter plus/minus, 1 for (1+x^2), -1 for (1-x^2)\n");
	scanf("%d",&pm);
	if(pm==1)
		sprintf(pmstr,"1+x^2");
	if(pm==-1)
		sprintf(pmstr,"1-x^2");
		
	fprintf(OutFile,"\t\t*****+++*****\n\n\n1pX2q.cpp\n\nusing x =");
	crstrXlr(X,10,istring);
	fprintf(OutFile,"%s",istring);		
	
	
		
	if(a*b< 0 && xlt1 ==1)
	 {
		if(pm == -1)
			X.sgn *=-1;
		Divlxr( X,one,&X,Ndp);	
		a=abs(a);
		b=abs(b);
		xlt1 =2;		 
	 }		
		
	smx =tIx=one;
	fct  =one;
	prd = q;
	smIx = X;
	Kx2= one;
	prodxlr(X,X,&x2,Ndp);
	fprintf(OutFile,"\n\n +-x^2=");
	crstrXlr(x2,10,istring);
	fprintf(OutFile,"%s",istring);
	fprintf(OutFile,"\n\nusing q =a/b \t a= %ld b= %ld",a,b);
	fprintf(OutFile,"\n\n expansion of (%s) ^(%ld/%ld) =",pmstr,a,b );
	x2.sgn =pm;
	k2 =x2;int i=1;
	if(xlt1 == 1)
	{	
		Kx2 = x2;
		while(i<=a)
		{			
			prodxlr(x2,Kx2,&Kx2,Ndp);
			i++;
		}
		NRoot(Kx2,b,&Kx2,Ndp);			//		X^2q
			
			
		Divlxr( x2,one,&x2,Ndp);
	}
		
	 i=1;
	// if(xlt1 ==0 or xlt1 == 2)
	while (smIx.exp -tIx.exp<Ndp)
	//while(i < 10)
	{
		prodxlr(k2,prd,&termx,Ndp);
		Divlxr( fct,termx,&termx,Ndp);
		addXlr(smx,termx,&smx,Ndp);
 
		prodxlr(X,termx,&tIx,Ndp);
		itoxR(2*i+1,&B);
		Divlxr( B,tIx,&tIx,Ndp);
		addXlr(smIx,tIx,&smIx,Ndp);
		if(i/100*100 ==i)
		{
			crstrXlr(termx,60,istring);
			printf("\nterm\n%s",istring);
			
			crstrXlr(tIx,60,istring);
			printf("\ntIx Int\n%s",istring);
			printf("\n\t%d",i);
		}
		i++;			
		itoxR(i,&bufi);	
		prodxlr(bufi,fct,&fct,Ndp);			 
		prodxlr(x2,k2,&k2,Ndp);
		
		itoxR(1-i,&bufi);
		addXlr(q,bufi,&bufi,Ndp);	
		prodxlr(bufi,prd,&prd,Ndp);
		/*if(i/10*10 ==i)
		{			
			crstrXlr(k2,60,istring);
			printf("\n i  %d\n k2\n%s",i,istring);
			
			crstrXlr(prd,Ndp,istring);
			printf("\nprd Int\n%s",istring);
			
		}*/
	}
	
	
	time(&t2);
	tsec=difftime(t2, t1);	
	
	crstrXlr(smx,Ndp,istring);
	fprintf(OutFile,"\n%s",istring);
	
	crstrXlr(smIx,Ndp,istring);
	fprintf(OutFile,"\n\nand\n\nintegral val \n%s",istring);
	
	itoxR(6,&bufi);					
	prodxlr(bufi,smIx,&bufi,Ndp);
	crstrXlr(bufi,Ndp,istring);
	fprintf(OutFile,"\nPI estimate\n%s",istring);
      
	fprintf(OutFile,"\r\n\n total compute time  %lf sec\n",tsec);

}
