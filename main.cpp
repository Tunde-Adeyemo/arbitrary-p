
/*
Place this file in the same directory as functns.cpp to compile.
main(), tests the routines in functnw.cpp and writes the output to local file tstfns.txt, it uses invtan() and
the Gregory series to efficiently calculate Pi to the desired precision then it uses NRoot() to produce root(1/2),
knowing cos(Pi/4) = root(1/2) it then calculates cos(Pi/4) - root(1/2) and writes the value. Finally it writes the 
number of seconds used for computation.
*/

#include "./functns.cpp"
int main(int argc, char* argv[])
{
	//FILE* OutFile;
 	NxlongR f,f2,r5,sin,cos,Pi4;
	char istring[2200];
	
	int Ndp=500; //  Ndp max ~ncells*6*0.6
	
	
	NxlongR one;	
	itoxR(1,&one);

	
	double tsec;

	time_t t1,t2;
        char filename[25];
        
           strcpy(filename, ".//tstfns.txt");


    FILE* OutFile = fopen(filename, (char *) &"w");
    if (OutFile == NULL) {
        printf("Cannot create output files\n");
        printf("Usage: Execute native/seq_demo2 from ..."); 
        exit(0);
    }
    
	time(&t1);
	fprintf(OutFile,"\n\n  calculate Pi =16 arctan 0.2-4 arctan(1/239)\n");


	itoxR(5,&r5);
	Divlxr( r5,one,&r5,Ndp);
	invtan(r5,&f,Ndp);
	itoxR(16,&r5);
	prodxlr(r5,f,&f,Ndp);

	itoxR(239,&r5);
	Divlxr( r5,one,&r5,Ndp);
	invtan(r5,&f2,Ndp);
	itoxR(4,&r5);
	r5.sgn=-1;
	prodxlr(r5,f2,&f2,Ndp);		
	addXlr(f,f2,&f,Ndp);		//	*** Pi


			
	crstrXlr(f,Ndp,istring);
	fprintf(OutFile,"\r\n\nPI =\r%s",istring);

	itoxR(4,&r5);
	Divlxr( r5,f,&r5,Ndp);// Pi/4
	Pi4=r5;	
	//CosSin(mpf_t x,int cs,mpf_t *CoS)	//cs 0,cos or 1,sin 
	CosSin(Pi4,0,Ndp,&cos);
	crstrXlr(cos,Ndp,istring);
	fprintf(OutFile,"\r\n\ncos(PI/4) =\r%s",istring);
	
	itoxR(2,&r5);
	Divlxr( r5,one,&r5,Ndp);
	NRoot(r5,2,&r5,Ndp);
	r5.sgn=-1;
	addXlr(r5,cos,&r5,Ndp);
	
	crstrXlr(r5,20,istring);
	fprintf(OutFile,"\r\n\ncos(PI/4) - rt(0.5) =\n%s",istring);
	
	time(&t2);
	tsec=difftime(t2, t1);
	fprintf(OutFile,"\r\n\n total compute time  %lf sec\n",tsec);    
	
}


