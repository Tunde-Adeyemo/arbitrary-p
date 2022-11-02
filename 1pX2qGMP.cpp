# include <stdlib.h>
#include <math.h>
#include <cstdio>
#include <string.h>
#include <time.h>
#include <gmp.h>



/*Developed by Olatunde Adeyemo �2022
 * from 1pX2q.cpp to obtain integral of fns of type (1+-x^2)^q

let us assume we are working to n dp through int Ndp
*/
//#define USEC_TO_SEC 1.0e-6

//using namespace std;
//namespace std does not work with GMP - forces error

//ifstream InFile;   // used with  rdXlr()

#include "./fns.cpp"
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
	int Ndp=5000;
	//FILE* OutFile;
 	//NxlongR one,X,x2,q,B,fct,prd,k2,termx,tIx,smx,smIx,bufi,x2q;
 	mpf_set_default_prec ( Ndp*log(10.0)/log(2.0));
 	mpf_t one, X,x2,q,B,fct,prd,k2,termx,tIx,smx,smIx,bufi,x2q,dm2,Kx2;
 	
 	mpf_init_set_d(one,1.0); 	
 	mpf_init(X);
 	mpf_init(x2);
 	mpf_init(q);
 	mpf_init(B);
 	mpf_init(fct);
 	mpf_init(prd);
 	mpf_init(k2);
 	mpf_init(termx);
 	mpf_init(tIx);
 	mpf_init(smx);
 	mpf_init(smIx);
 	mpf_init(bufi);
 	mpf_init(x2q);
 	mpf_init(dm2);
 	mpf_init(Kx2);
 	
	char istring[8000], pmstr[10];
	
	int  xlt1=0,pm;
	long a,b,expo;
	double Xdbl;
	
	double tsec;

	time_t t1,t2;
        char filename[25]; 	
    
     strcpy(filename, ".//1pX2qG.txt");


    FILE* OutFile = fopen(filename, (char *) &"w");
    if (OutFile == NULL) {
        printf("Cannot create output files\n");
        printf("Usage: Execute native/seq_demo2 from ..."); 
        exit(0);
    }
    setbuf(OutFile,NULL );
    
	time(&t1);
	
	
	//itoxR(1,&one);
    
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
	scanf("%li",&expo);
	//xp= pow(10,expo);   //// wrong accuracy !!!ONLY DOUBLE
	
	pm=5;	
	if (expo < 0){
	 pm =-5;
	 expo *= -1;	
	}
	unsigned long xp = expo; 
	itoxR(10,&B);
	mpf_pow_ui (dm2,B,xp);
		
	
	mpf_set_d( x2,Xdbl);	
	if(pm == 5)
		mpf_mul (X,x2,dm2);
	else
		mpf_div (X,x2,dm2);
	
	mpf_get_d_2exp (&b,X);
	//printf(" exp X %ld",b);	
	if (b > 0)	// 0.1 and decimals give exp 0
		xlt1=1;
	printf("\nenter rational value q  in the form of 2 integer entries a and b\n \
	 q = a/b  b may be 1,\nenter a\n");
	 scanf("%ld",&a);
	 printf("enter b\n");
	 scanf("%ld",&b);
	 itoxR(a,&bufi);
	 itoxR(b,&B);
	 Divlxr( B,bufi,&q);        // ** "real" q used in (1-u)^q	
	 
	 
	printf("\nenter plus/minus, 1 for (1+x^2), -1 for (1-x^2)\n");
	scanf("%d",&pm);
	if(pm==1)
		sprintf(pmstr,"1+x^2");
	if(pm==-1)
		sprintf(pmstr,"1-x^2");
		
	fprintf(OutFile,"\t\t*****+++*****\n\n\n1pX2qGMP.cpp\n\nusing x =");
	crstrXlr(X,10,istring);
	fprintf(OutFile,"%s",istring);	
	
	
		
		
	
	mpf_set(smx,one);
	mpf_set(tIx,one);
	mpf_set(fct,one);
	mpf_set(Kx2,one);
	
	/*printf("X ");
	mpf_out_str (stdout, 10,20,X);
	printf("\n q \t");
	mpf_out_str (stdout, 10,25,q);
			printf("\n ");*/
	
	mpf_set(prd,q);
	
	mpf_set(smIx,X);	
	prodxlr(smIx,X,&x2);
	fprintf(OutFile,"\n\n (+-)x^2= ");	
	crstrXlr(x2,10,istring);
	fprintf(OutFile,"%s ",istring);
	fprintf(OutFile,"\n\nusing q =a/b \t a= %ld b= %ld",a,b);
	fprintf(OutFile,"\n\n expansion of (%s) ^(%ld/%ld) =",pmstr,a,b );
	int i=1;
	if(xlt1 == 1)
	{	
		mpf_set(Kx2,x2);
		while(i<=a)
		{
			mpf_set(bufi,Kx2);
			prodxlr(x2,bufi,&Kx2);
			i++;
		}
		mpf_set(bufi,Kx2);
		NRoot(bufi,b,&Kx2);			//		X^2q
		if (pm == -1)
			mpf_neg(Kx2,Kx2);	
			
		mpf_set(bufi,x2);
		Divlxr( bufi,one,&x2);
	}
	
	if (pm == -1)
		mpf_neg(x2,x2);	
	
	mpf_set(k2,x2);	
			
	
	 	
	long exsmI,extI;
	mpf_get_d_2exp (&exsmI,smIx);	
	mpf_get_d_2exp (&extI,tIx);
	
	while (exsmI -extI<Ndp*log(10.0)/log(2.0))
	//while(i < 10)
	{			
		prodxlr(k2,prd,&B);					
			
		Divlxr( fct,B,&termx);
		mpf_set(B,smx);
		addXlr(B,termx,&smx);
 
		prodxlr(X,termx,&tIx);
		itoxR(2*i+1,&B);
		mpf_set(bufi,tIx);
		Divlxr( B,bufi,&tIx);
		addXlr(smIx,tIx,&smIx);
		if(i/100*100 ==i)			//
		{
			crstrXlr(termx,60,istring);
			printf("\t termx \n%s\n",istring);
			//mpf_out_str (stdout, 10,30,termx);
			crstrXlr(tIx,60,istring);
			printf("\ntIx Int\n%s\n",istring);
			//mpf_out_str (stdout, 10,30,tIx);
			printf ("\n exsmI -extI %ld \n",exsmI -extI);
			printf("\n\t%d",i);
		}//
		i++;			
		itoxR(i,&bufi);	
		mpf_set(B, fct);
		prodxlr(bufi,B,&fct);	
		mpf_set(B,k2);		 
		prodxlr(x2,B,&k2);
		
		itoxR(1-i,&B);
		addXlr(q,B,&bufi);	
		mpf_set(B, prd);
		prodxlr(bufi,B,&prd);
		/*if(i/10*10 ==i)
		{			
			crstrXlr(k2,60,istring);
			printf("\n i  %d\n k2\n%s",i,istring);
			
			crstrXlr(prd,Ndp,istring);
			printf("\nprd Int\n%s",istring);
			
		}*///
		mpf_get_d_2exp (&exsmI,smIx);	
		mpf_get_d_2exp (&extI,tIx);
	}
	mpf_mul(smx,Kx2,smx);
	mpf_mul(smIx,Kx2,smIx);
	time(&t2);
	tsec=difftime(t2, t1);
	//crstrXlr(smx,Ndp,istring);
	//fprintf(OutFile,"\nfn of(1+x)\n%s",istring);
	/*prodxlr(smx,smx,&B,Ndp);
	crstrXlr(B,Ndp,istring);
	fprintf(OutFile,"\n(1+x)^2\n%s\n",istring);*
	*/
	
	crstrXlr(smx,Ndp,istring);
	fprintf(OutFile,"\n%s",istring);
	
	crstrXlr(smIx,Ndp,istring);
	fprintf(OutFile,"\n\nand\n\nintegral val \n%s",istring);
	
	itoxR(6,&bufi);			
	mpf_mul(B,bufi,smIx);
	crstrXlr(B,Ndp,istring);
	fprintf(OutFile,"\n6*sum int \n%s",istring);
      
	fprintf(OutFile,"\r\n\n total compute time  %lf sec\n",tsec);

}
