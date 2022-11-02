/* 
 * mar 21
 * Draft Newton method for PI using expansion of (1+x^2)^0.5
 * based on Newtroot.cpp
 * improved from Newtpigmp.cpp
 * 
 */	

	
	# include <stdlib.h>
	# include <math.h>
	# include <stdio.h>
	# include <string.h>	
	# include <time.h>
	# include <gmp.h>	



inline void crstrXlr(mpf_t xlr,int sf,char* ostr)
	{
		char var[30], flstr[10000];
		int i,ph=0,i2;
		long int xp,ix;			
		
		ix = mpf_get_si(xlr);		
		if (mpf_cmp_si(xlr,ix) == 0)
		{			
			sprintf(ostr,"\r  %ld",ix);
			return;
		}
		mpf_get_str ( flstr,&xp, 10,sf+10,xlr);		
		
		if (flstr[0] == '-')
		{
			strncpy(var,flstr,2);
			var[2] = '.';
			ph =2;
		}
		else
		{
			strncpy(var,flstr,1);
			var[1] = '.';
			ph =1;
		}
		i2=ph;
		for(i=i2; i<i2+5;i++)
		{
			var[ph+1]=flstr[ph];
			ph++;
		}
		
		sprintf(ostr,"\r  %s",var);
		
		
			i2 =ph;
			for(i=2;i<=(sf+4)/6;i++)
			{	
				strncpy(var,&flstr[ph],6);
				strcat(ostr," ");
				strncat(ostr,var,6);		//strncat to accurately only copy 6 chars
				
				if(i/10*10==i)
					strcat(ostr,"\n  ");
				ph +=6;
			}
			sprintf(var," exp %li",xp-1);
			strcat(ostr,var);

	}


int main(int argc, char* argv[])
{	
	mpf_set_default_prec ( 50000);		// needs to be here to work
	int i;
	long int px,tx,dp;
	//mpf_t bincf[itr], nfonmbf[itr],x2,x,term[itr],Pi;
	mpf_t bincf, nfonmbf,x2,x,term,Pi, term0, bf;
	FILE* OutFile;
	char filename[25]	,Os[20000];
	strcpy(filename, ".//NewtPi.txt");
	OutFile = fopen( filename,"w");	
	
	long tsec;   
	time_t t1,t2;
		
	mpf_init(bincf);
	mpf_init(nfonmbf);
	mpf_init(term);
	mpf_init (term0);
	
	
	mpf_init_set_d(x2,-0.25);
	mpf_init_set_d(x,0.5);
	mpf_init_set_d(Pi,0.5);
	mpf_init_set_d(bf,1.0);	
	
	
	mpf_set_d( bincf,0.25);			// get area of triangle first sub from storage Pi
	mpf_set_d( term,0.75);				// we can now use them as dummy storage vars
	mpf_sqrt( term,term);
	mpf_mul( term,bincf,term);
	mpf_sub( Pi, Pi, term),	
	
	mpf_set_ui( bincf,1);
	mpf_div_ui(term0,bincf,2);
	mpf_set( bincf,term0);
	mpf_set( nfonmbf,term0);
	
	//term[1]=(x^3)*nfonmbf[1]/3; //integrate first term	
	mpf_mul(x,x2,x);
	mpf_set( term,term0);
	mpf_div_ui(term,term,3);
	//mpf_mul(term,term,nfonmbf);
	mpf_mul(term,term,x);	
	mpf_add(Pi,Pi,term);
	mpf_get_str ( Os,&px, 2,10,Pi);
	mpf_get_str ( Os,&tx, 2,10,term);
	i=2;
	dp = mpf_get_default_prec();
	time(&t1);
	
	while(px-tx < dp)
	{
		//nfonmbf[i]=nfonmbf[i-1]*(1.5-i);
		mpf_t dum;
		mpf_init_set_d( dum, 1.5-i);
		mpf_mul(nfonmbf,nfonmbf,dum);
		
		//bf *= i;
		mpf_mul_ui(bf,bf,i);
		mpf_div(bincf,nfonmbf,bf);
		mpf_mul(x,x,x2);
		
		mpf_mul( term,bincf,x);
		mpf_div_ui(term,term,2*i+1);				
		mpf_add(Pi,Pi,term);
		mpf_get_str ( Os,&tx, 2,10,term);
		i++;		
		
		/*if(i > itr-6)
		{
			//OutFile <<i<<" \t"<<nfonmbf[i]<<" \t"<<bincf[i]<<" \t"<<term[i]<<" \t\n";
			fprintf(OutFile,"\r\n %d \t",i);
			mpf_out_str(OutFile,10,10,nfonmbf);
			
			fprintf(OutFile," \t ");
			mpf_out_str(OutFile,10,10,bincf);
			
			fprintf(OutFile," \t ");
			mpf_out_str(OutFile,10,10,term);
		}	*/
	}
	time(&t2);
	tsec=difftime(t2, t1);
	printf("\n time section 1 is %ld s\n",tsec);
	
	printf("\r\n %d \t",i);
	mpf_out_str(stdout,10,10,nfonmbf);
			
	printf(" \tbincf ");
	mpf_out_str(stdout,10,10,bincf);
	
	printf(" \tterm ");
	mpf_out_str(stdout,10,10,term);
	
		
	printf("\n\t");
	
	mpf_mul_ui( Pi, Pi, 12);	
	mpf_out_str(stdout,10,100,Pi);
	fprintf(OutFile,"\n");
	
	mpf_out_str(OutFile,10,20,Pi);
	fprintf(OutFile,"\n\t last term is \n\t");
	mpf_out_str(OutFile,10,20,term);
	
	crstrXlr(Pi,600,Os);
    fprintf(OutFile,"\r\n\n\testimate for PI\r \t%s",Os);
}
	
