/*
 *
 *
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
//using namespace std;*/
   
    

	// addXlr 06/01/14
	void addXlr(mpf_t a, mpf_t b,mpf_t* c)
	{		
		mpf_add( *c,a,b);
	}
	
	inline void prodxlr(mpf_t a, mpf_t b,mpf_t* c)
	{
		mpf_mul(*c,a,b);

	}
	
	void Divlxr(mpf_t divisor,mpf_t divend,mpf_t* quot)
	{
		mpf_div(*quot,divend,divisor);	
	}
	
	void itoxR(int intg, mpf_t* xlo)
	{				
		mpf_init_set_si(*xlo,intg);
	}

	void invtan(mpf_t x,mpf_t* Itan)
	{
		mpf_t xp,div,term,two ;
		
		mpf_init( xp);
		mpf_init( div);
		mpf_init( term);
		mpf_init( two);
		long dp= mpf_get_default_prec(),xxp, txp;
		
		int sg=1, sdp=-dp;
		
		mpf_get_d_2exp (&xxp,x);
			if(xxp<0)
				sdp=-dp+xxp;

				
		mpf_set (xp,x);
		mpf_set (term,x);
		mpf_set (*Itan,x);
		itoxR(2,&two);
		itoxR(1,&div);
		mpf_get_d_2exp (&txp,term);		
		mpf_mul(x,x,x);
		
		while(txp>sdp)
		{
			sg*=-1;			
			mpf_add(div,div,two);
			prodxlr(xp, x,&xp);
			Divlxr( div,xp,&term);
			
			if (sg>=0)				
				mpf_add(*Itan,*Itan,term);			
			else
				mpf_sub(*Itan,*Itan,term);
			
			mpf_get_d_2exp (&txp,term);	
		}
		return;

	}
	void CosSin(mpf_t x,int cs,mpf_t *CoS)	//cs 0,cos or 1,sin 
	{
		mpf_t xn,fac,f2,t;
		mpf_init(fac);
		mpf_init(f2);
		mpf_init(xn);
		mpf_init(t);
		long fct,sg;
		long Ndp = mpf_get_default_prec(),xf2;
		
		itoxR(1,CoS);
		
		if(cs==1)	
		mpf_set(*CoS,x);		
	
		mpf_set(xn,*CoS);

		mpf_mul(x,x,x);
		fct=2;
		if(cs==1)
		fct=3;

		itoxR(2,&fac);
		if(cs==1)
		itoxR(6,&fac);
		
		//mpf_get_d_2exp (&xf2,f2);
		xf2 =0;
		sg=- mpf_sgn (*CoS); // when x -ve		
		while (xf2 >-Ndp)
		{
				
			mpf_mul(xn,x,xn);
			
			mpf_div(f2,xn,fac);
			
			mpf_abs ( f2, f2 );
			if (sg <0)
			mpf_neg(f2,f2);
			
			for(int i=1;i<=2;i++)
			{
				fct++;
				itoxR(fct,&t);
				
				mpf_mul(fac,fac,t);
			}		

			
			mpf_add(*CoS,*CoS,f2);
			sg*=-1;
			
			mpf_get_d_2exp (&xf2,f2);
						
		}		
		return;
	}


	inline void crstrXlr(mpf_t xlr,int sf,char* ostr)
	{
		char var[30], flstr[100000];      // flstr must match or greater than *ostr
		
		int i,ph=0,i2;
		long int xp,ix;		
		
		ix = mpf_get_si(xlr);		
		if (mpf_cmp_si(xlr,ix) == 0)
		{			
			sprintf(ostr,"\n  %ld",ix);
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
		for(i=i2; i<=i2+5;i++)				// *****
		{
			var[ph+1]=flstr[ph];
			ph++;
		}
		
		sprintf(ostr,"\n   ");
		strncat(ostr,var,ph);
				
			i2 =ph;      
			for(i=2;i<=(sf+4)/6;i++)
			{	
				strncpy(var,&flstr[ph-1],6);		//********
				strcat(ostr," ");
				strncat(ostr,var,6);		//strncat to accurately only copy 6 chars
				
				if(i/10*10==i)
					strcat(ostr,"\n  ");
				ph +=6;
				if(i/100*100==i)
				{
					sprintf(var, " %i\n  ",(i*6)/100);
					strcat(ostr,var) ;  
				}
			}
			sprintf(var," exp %li",xp-1);
			strcat(ostr,var);

	}



	void NRoot(mpf_t A, int N,mpf_t* xf )		//,int prc)
	{
		
	double x=1.0,rm,aA = mpf_get_d(A);
	int n=2,i;
	long int mi,aexp ;
	
	long int prc =  mpf_get_default_prec (  );
	
	mpf_t ff,ew;
	
//		Set the default precision to be at least prec bits. All
	
	mpf_init_set_d ( *xf, x);
			
	mpf_init (ew);
	mpf_init (ff);
	rm = fabs(pow(x,n)-aA);
	mpf_set_d(ew,rm);
		
	mpf_get_d_2exp (&mi,ew);
	mpf_get_d_2exp (&aexp,A);
	while(mi >-prc+aexp)
		{
				
			//  a(i+1) = a(i)/(1+(a(i)^n - A)/nA)		
			mpf_set(ff,*xf);
			mpf_set(ew,A);		
			for (i=1;i<n;i++)	
				mpf_mul( ff, *xf,ff );
			
			mpf_sub(ew,ff,ew);	// ew ,y
			
			mpf_set(ff,A);
			mpf_mul_ui(ff,ff,n);		
			mpf_div(ew,ew,ff);
			
			mpf_add_ui(ew,ew,1);
			mpf_div(*xf,*xf,ew);
			
			mpf_set(ff,*xf);
			mpf_set(ew,A);		
			for (i=1;i<n;i++)	
				mpf_mul( ff, *xf,ff );
			
			mpf_sub(ew,ff,ew);		
			mpf_get_d_2exp (&mi,ew);			
					
		}
	
	}




