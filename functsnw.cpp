	
	# include <stdlib.h>
	# include <math.h>
	# include <stdio.h>
	# include <string.h>
	# include <iostream>	
	# include <time.h>
	# define ncells 4000 //remember you need extra cells for accuracy in division using Divlx

/*Developed by Olatunde Adeyemo �2013
PIsf3 is derived from expo.cpp to calculate exp(x)

let us assume we are working to n dp through int Ndp
*/
//#define USEC_TO_SEC 1.0e-6

using namespace std;


/*double wallTime() {
    double seconds;
    struct timeval tv;

    gettimeofday(&tv, NULL);

    seconds = (double) tv.tv_sec; // seconds since Jan. 1, 1970
    seconds += (double) tv.tv_usec * USEC_TO_SEC; // and microseconds
    return seconds;
}*/

struct NxlongR {
	int mant[ncells]; //each cell of array mant carries 6 digits of our extended real number
	int exp; //exp is the exponent of our first non zero digit of mant[1]
	int sgn; //overall sign
};


/*			may 22
A compilation of the most updated functions of ArithXlr
includes corrected addXlr() , Divlxr(),improved NRoot().
Basic functions  CosSin(), invtan() and itoxR() are included
also the print function crstrXlr() in this form it displayed dp/100 in the line between blocks
rdXlr(), can read an NxlongR variable written by crstrXlr from a file.

inline void addXlr(NxlongR a, NxlongR b, NxlongR * c, int dp)
inline void prodxlr(NxlongR a, NxlongR b,NxlongR* c,int dp) 
inline void Divlxr(NxlongR divisor, NxlongR divend, NxlongR * quot, int dp)
void itoxR(int intg,NxlongR* xlo)
void invtan(NxlongR x,NxlongR* Itan,int dp)
inline void crstrXlr(NxlongR xlr,int sf,char* ostr)
void CosSin(NxlongR x,int cs,int Ndp,NxlongR* CoS)	//cs 0,cos or 1,sin 
void NRoot(NxlongR A, int N,NxlongR* x,int Ndp) // revised improved 28/03/21  remodelled 05/22
void rdXlr(char fn[25],NxlongR* x,int mntit,int ofl)
*/

//int df=0; 		//diagnostic global variable


void rippleup(NxlongR* a)
{
	for(int i=ncells-1;i>=2;i--)
		if(a->mant[i]>=1e6)
		{
			int carry=a->mant[i]/1e6;
			a->mant[i-1]+=carry;
			a->mant[i]-=carry*1e6;
		}
}

inline void xnegc(NxlongR * c, int dp) {
	double Dbl;
	int i, i2, ep, shft, shft2;
	for (i = ncells - 1; i >= 2; i--) // check for carry-up and neg results
		if (c -> mant[i] < 0) {
			c -> mant[i - 1] --;
			c -> mant[i] += 1000000;
		}
	for (i = ncells - 1; i > 1; i--) {
		if (c -> mant[i] >= 1e6)
			rippleup(c);
	}

	if (c -> mant[1] >= 1e6) {
		Dbl = (double) c -> mant[1];
		ep = 0;

		while (Dbl >= 1e6) {
			Dbl /= 10;
			ep++;
		}
		c -> mant[1] = Dbl;
		Dbl -= (double) c -> mant[1];

		for (i = 2; i <= (c -> exp + dp) / 6 + 1; i++) {
			int j;
			for (j = 1; j <= 6 + ep; j++)
				Dbl *= 10;

			Dbl += c -> mant[i];
			for (j = 1; j <= ep; j++)
				Dbl /= 10;

			c -> mant[i] = Dbl;
			Dbl -= (double) c -> mant[i];
		}
		c -> exp += ep;
	}

	ep = 0;
	i = 1;
	if (c -> mant[1] == 0) {
		i = 2;
		while (c -> mant[i] == 0)
			i++;
		for (i2 = 1; i2 < dp/ 6 + 2-i; i2++) {
			c -> mant[0] = c -> mant[i + i2 - 1];
			c -> mant[i2] = c -> mant[0];
		}


		if (i * 6 >= dp) {

			c -> sgn = 0;
			return;
		}
	}

	while (c -> mant[1] < 1e5) {
		ep++;
		c -> mant[1] *= 10;
	}	
	c -> exp -= ep + (i - 1) * 6;
	if (ep > 0) {


		shft2 = shft = c -> mant[2];

		for (i2 = 1; i2 <= 6 - ep; i2++)
			shft /= 10;

		c -> mant[1] += shft;


		for (i2 = 1; i2 <= 6 - ep; i2++)
			shft *= 10;

		shft = shft2 - shft;


		i2 = 2;
		while (i2 * 6 <= dp) {
			int i3;

			for (i3 = 1; i3 <= ep; i3++)
				shft *= 10;
			c -> mant[i2] = shft;


			shft = shft2 = c -> mant[i + i2];

			for (i3 = 1; i3 <= 6 - ep; i3++)
				shft /= 10;
			c -> mant[i2] += shft;

			for (i3 = 1; i3 <= 6 - ep; i3++)
				shft *= 10;
			shft = shft2 - shft;

			i2++;

		}
	}

}

// addXlr 06/01/14 // 25/04/21 correction for a ~= b and mant[2]
inline void addXlr(NxlongR a, NxlongR b, NxlongR * c, int dp) {
	NxlongR b2;	
	int bign = a.exp - b.exp, xpb6 = bign / 6, i, i2, shft,i2m;
	//double Dbl;

	if (a.sgn == 0) { * c = b;
		return;
	}
	if (b.sgn == 0) { * c = a;
		return;
	}

	for (i = 1; i < ncells; i++)
		c -> mant[i] = 0;


	if (bign > 0) {
		for (i = 1; i <= xpb6; i++)
			c -> mant[i] = a.mant[i];

		int rem = bign - xpb6 * 6;
		b2.mant[1] = 0;
		for (i = 1; i <= dp / 6 + 1; i++) // b must be aligned before adding
		{
			shft = b.mant[i];
			for (i2 = 1; i2 <= rem; i2++)
				shft /= 10;

			b2.mant[i] += shft;
			for (i2 = 1; i2 <= rem; i2++)
				shft *= 10;

			shft = b.mant[i] - shft;
			for (i2 = 1; i2 <= 6 - rem; i2++)
				shft *= 10;
			b2.mant[i + 1] = shft;

		}
		for (i = dp / 6 + 1; i > xpb6; i--)


			switch (a.sgn * b.sgn) {
			case 1:
				{
					c -> mant[i] = a.mant[i] + b2.mant[i - xpb6];
					break;
				}

			default:
				{
					c -> mant[i] = a.mant[i] - b2.mant[i - xpb6];
					break;
				}


		}
		c -> exp = a.exp;
		c -> sgn = a.sgn;
	}

	if (bign < 0) {
		bign = -bign;
		xpb6 = bign / 6;

		for (i = 1; i <= xpb6; i++)
			c -> mant[i] = b.mant[i];
		int rem = bign - xpb6 * 6;

		b2.mant[1] = 0;
		for (i = 1; i <= +dp / 6 + 1; i++) // a must be aligned before adding
		{

			shft = a.mant[i];
			for (i2 = 1; i2 <= rem; i2++)
				shft /= 10;

			b2.mant[i] += shft;
			for (i2 = 1; i2 <= rem; i2++)
				shft *= 10;

			shft = a.mant[i] - shft;
			for (i2 = 1; i2 <= 6 - rem; i2++)
				shft *= 10;
			b2.mant[i + 1] = shft;


		}
		for (i = dp / 6 + 1; i > xpb6; i--)


			switch (a.sgn * b.sgn) {
			case 1:
				{
					c -> mant[i] = b.mant[i] + b2.mant[i - xpb6];
					break;
				}

			default:
				{
					c -> mant[i] = b.mant[i] - b2.mant[i - xpb6];
					break;
				}
		}

		c -> exp = b.exp;
		c -> sgn = b.sgn;


	}


	if (bign == 0) {
		i = 0;
		if (a.sgn * b.sgn == 1){
			while (i * 6 <= dp) 
                        {
				c -> mant[i] = b.mant[i] + a.mant[i];
				++i;
                        }

			}
                c -> sgn = b.sgn;
                c -> exp = b.exp;
                
		if (a.sgn * b.sgn == -1) {

			while (a.mant[i + 1] == b.mant[i + 1])
				i++;
			
			if (i * 6 >= dp) {
				c -> mant[1] = 0;
				c -> sgn = 0;
				return;
			}
                        
                
						
			int i2 = i2m=(dp+5)/6;
			while(i2 >=1)

				{
                                    if (a.mant[i + 1] > b.mant[i + 1])
                                            c -> mant[i2] = a.mant[i + i2] - b.mant[i + i2];						
                                    else					
                                            c -> mant[i2] = b.mant[i + i2] - a.mant[i + i2];

											
					if (c-> mant[i2+1] < 0 && i2+1 <=i2m)
					{
						c-> mant[i2] --;
						c-> mant[i2+1] += 1000000;						
					}	

					i2 --;					
				}
				if (a.mant[i + 1] > b.mant[i + 1])
					{						
						c -> exp = a.exp - i * 6;
						c -> sgn = a.sgn;
					}
					else
					{						
						c -> exp = b.exp - i * 6;
						c -> sgn = b.sgn;
					}						
			}
	}
        int im =2;
        if ( c->mant[1] == 0)
	{
		while(c->mant[im] == 0)
			im++;
		if ( c->mant[im] < 0)
		{	
			c->sgn *= -1;
			for(int im2 = i2m;im2 >= im; im2--)
                        {
                            c ->mant[im2] *= -1;

                            if (c-> mant[im2+1] < 0 && im2 < i2m)
                                            {
                                                    c-> mant[im2] --;
                                                    c-> mant[i2+1] += 1000000;						
                                            }
                        }
			
		}
	}
	xnegc(c, dp);
	
	if (c -> mant[1] == 0)
		c -> sgn = 0;
}

inline void prodxlr(NxlongR a, NxlongR b,NxlongR* c,int dp) //25/03/14
{
	int i,j,sf1b6,sf2b6,f,m,bk,ep=0;

	double Dbl,man1;
	sf1b6=(dp+5)/6+1;	//(a.exp+dp)/6+1;
	sf2b6=(dp+5)/6+1;	//(b.exp+dp)/6+1;
	if(sf1b6<1)
		sf1b6=1;
	if(sf2b6<1)
		sf2b6=1;

	man1=(double)a.mant[1]*b.mant[1]/1e10;

	while(man1>=10)
	{
		man1/=10;
		ep++;
	}
	c->sgn=a.sgn*b.sgn;
	c->exp=a.exp+b.exp+ep;  //need to determine product of man[1] for magnitude
	for(i=1;i<ncells;i++)
	c->mant[i]=0;
	
	if(a.sgn==0)
	{
		c->sgn=0;
		return;
	}

	if(b.sgn==0)
	{
		c->sgn=0;
		return;
	}

	for(i=1;i<sf1b6;i++)           //for(i=1;i<=sf1b6;i++)
	for(j=1;j<sf2b6;j++)
	{
		Dbl=0.0;
		if ((i+j+1)*6<dp)
			Dbl=(double)a.mant[i]*b.mant[j];
		f=0;
		m=0;
		bk=0;

		if(Dbl>=1e11)
		{
			f=Dbl/1e11;
			Dbl-=f*1e11;
		}
	
		if(Dbl>=1e5)
		{
			m=Dbl/1e5;
			Dbl-=m*1e5;
		}
	
		bk=Dbl*10;
	
		if(i==1&&j==1)
		switch(ep)
		{

			case 1:
			{
				c->mant[1]=f;
				c->mant[2]=m;
				c->mant[3]=bk;
				break;
			}

			default:
			{
				//c.mant[1]=f;
				c->mant[1]=m;
				c->mant[2]=bk;
				break;
			}
		}
	
		
		if((i!=1||j!=1)&&((i+j)*6<=dp))
		{


			switch(ep)
			{

				case 1:
				{
					c->mant[i+j-1]+=f;
					c->mant[i+j]+=m;
					c->mant[i+j+1]+=bk;
					break;
				}

				default:
				{
					c->mant[i+j-2]+=f;
					c->mant[i+j-1]+=m;
					c->mant[i+j]+=bk;
					break;
				}
			}
					
		}
	}
	if(c->mant[1]==0)
	{
		c->sgn=0;
		return;
	}

	
		if (ep>0)
		{
			int ex=0,shft,shft2,i3;
			while(c->mant[1]<1e5)
			{
				ex++;
				c->mant[1]*=10;
			}
							
			shft2=shft=c->mant[2];
                        int i2;

			for( i2=1;i2<=6-ex;i2++)
				shft/=10;

			c->mant[1]+=shft;		

		
			for(i2=1;i2<=6-ex;i2++)
				shft*=10;

			shft=shft2-shft;


			i2=2;
			while(i2*6<=dp)
			{

				for( i3=1;i3<=ex;i3++)
					shft*=10;
				c->mant[i2]=shft;
				

				shft=shft2=c->mant[1+i2];			

				for( i3=1;i3<=6-ex;i3++)
				shft/=10;
				c->mant[i2]+=shft;

				for( i3=1;i3<=6-ex;i3++)
				shft*=10;
				shft=shft2-shft;

				i2++;

			}
		}
	
	
	xnegc(c,dp);

}



//09/01/14      20/03/16  // oscill problem  solved in rtNA.cpp in ./workbook
inline	void Divlxr(NxlongR divisor,NxlongR divend,NxlongR* quot,int dp)
	{

		NxlongR rem,Int2,quot2,lpb,one ;
		//char str[100];
		double Dvisor,Dvend,Dquot;
		
		int ex, rmt,rp =1 ;
		dp*=1.2;
				
		if(divisor.sgn==0)
		{
			quot->exp=12345;
			quot->sgn=divisor.sgn*divend.sgn;
			return;
		}
		if(divend.sgn==0)
		{
			*quot= divend;
			quot->sgn=0;
			return;
		}


		Dvisor= divisor.mant[1]*1e-5+divisor.mant[2]*1e-11;
		Dvend= divend.mant[1]*1e-5+divend.mant[2]*1e-11;
		Dquot=Dvend/Dvisor;
		ex=0;
		while (Dquot>=10)
		{
			Dquot/=10;
			ex++;
		}
		while (Dquot<1.0)
		{
			Dquot*=10;
			ex--;
		}
		quot->exp=ex+divend.exp-divisor.exp;
		quot->sgn=divend.sgn*divisor.sgn;


		Dquot*=1e5;
		quot->mant[1]=Dquot;


		int i=1,sf=12;		//sf for quotient is different from remainder may require increasing for validity


		while((i+1)<=(sf+5)/6)
		{
			Dquot-=quot->mant[i];
			i++;
			Dquot*=1e6;
			quot->mant[i]=Dquot;
		}
		while(i<ncells-1)
		{
			i++;
			quot->mant[i]=0;

		}

		xnegc(quot,dp);
		
		//
		Int2 = divisor;   // divisor ~ divend
		Int2.sgn = -divend.sgn;
		addXlr(Int2, divend, &Int2, dp);
		Int2.sgn *=divend.sgn;
		
		if(Int2.exp - divend.exp <= -12)   // if divisor ~ divend  
		{  
			int tmp1;
			tmp1=1;
			
			for(i=1;i<ncells;i++)
				one.mant[i]=0;

			one.exp=0;
			one.sgn= 1;

			tmp1*=1e5;
			one.mant[1]=tmp1+0.05;
					
		    int qs = quot -> sgn;
		    
		    Dvend = Int2.mant[1] * 1e-5 + Int2.mant[2] * 1e-11 + Int2.mant[3] * 1e-17 ;
		    Int2.exp = Int2.exp - divend.exp;
		    
		    Dquot = Dvend / Dvisor; 
		    ex = 0;
		    while (abs(Dquot) >= 10) {
		            Dquot /= 10;
		            ex++;
		    }
		    while (abs(Dquot) < 1.0) {
		            Dquot *= 10;
		            ex--;
		    }
		    rem.exp = Int2.exp;
		    
		    Dquot *= 1e5;
		    rem. mant[1] = Dquot;

		    i = 1;
		    while ((i + 1) <= (sf + 5) / 6) {
		    Dquot -= rem. mant[i];
		    i++;
		    Dquot *= 1e6;
		    rem. mant[i] = Dquot;
		    }
		    while (i < ncells - 1) {
		            i++;
		            rem.mant[i] = 0; 
		    }

		    addXlr(one,  rem, quot, dp);		
		    quot->sgn = qs;
		} 			//   **

		prodxlr(*quot,divisor,&Int2,dp);
	            
        	Int2.sgn*=-1;
		addXlr(divend,Int2,&rem,dp);

		if(rem.sgn==0)
			return;

		rmt=rem.exp;
		rp=1;
	
		//int ii=1;    // for diag print
		while(divend.exp-rem.exp<0.82*dp)		// 0.82  dp =1.5 dp(in)0.68, 1.3 0.8
		{
	
			Dvend= rem.mant[1]*1e-5+rem.mant[2]*1e-11;
			Dquot=Dvend/Dvisor;
			ex=0;

				while (Dquot>=10)
				{
					Dquot/=10;
					ex++;
				}
				while (Dquot<1.0)
				{
					Dquot*=10;
					ex--;
				}


				quot2.exp=ex+rem.exp-divisor.exp;
				quot2.sgn=rem.sgn*divisor.sgn;

				Dquot*=1e5;
				quot2.mant[1]=Dquot;


				i=1;
				while((i+1)<=(sf+5)/6)
				{
					Dquot-=quot2.mant[i];
					i++;
					Dquot*=1e6;
					quot2.mant[i]=Dquot;
				}
				while(i<ncells-2)
				{
					i++;
					quot2.mant[i]=0;

				}
				xnegc(&quot2,dp);
			
				addXlr(*quot,quot2,quot,dp);

				prodxlr(*quot,divisor,&Int2,dp);
				Int2.sgn*=-1;
				addXlr(divend,Int2,&rem,dp);

				
				if(rem.sgn==0)
					return;
				
		
				if (rem.exp ==rmt)
				{
					rp ++;
					if (rp == 8)
						lpb = quot2;
										
					if(rp > 8 )
					{
	
						addXlr(lpb,quot2,&quot2,dp);
						if (quot2.sgn == 0)
								return;
						addXlr(*quot,quot2,quot,dp);
	
					}						
					
				}
				else
				{
					rmt = rem.exp;
					rp = 1;
				}
		}   // end while

	
	}

void itoxR(int intg,NxlongR* xlo)
{
	double tmp1;
	int ep=0,i;
	
	tmp1=abs(intg);
	xlo->sgn= intg/tmp1;	//tmp1 is already pos
	while(tmp1>=10.0)
	{
		tmp1/=10;
		ep++;
	}
	for(i=1;i<ncells;i++)
		xlo->mant[i]=0;

    xlo->exp=ep;	

	tmp1*=1e5;
	xlo->mant[1]=tmp1+0.05;
	i=1;
	int sf=14;
	tmp1-=(double)xlo->mant[1];
	while((i+1)*6<sf)
	{
		i++;
		tmp1*=1e6;
		xlo->mant[i]=tmp1;
		if((i+2)*6>sf)
		xlo->mant[i]=tmp1+0.05;// careful not to round down
		tmp1-=xlo->mant[i];
	}

}


inline void crstrXlr(NxlongR xlr,int sf,char* ostr)
{
	char var[30];
	int i,tprt,t2,i0,ep=0;
	tprt= t2=xlr.sgn*xlr.mant[1];
	
	for(i=1;i<=5;i++)
		tprt/=10;
	
	sprintf(ostr,"\n  %i.",tprt);
	
	//sprintf(ostr,"\n   ");
	//	strcat(ostr,tprt);
		for(i=1;i<=5;i++)
		tprt*=10;

		tprt=abs(t2-tprt);

			while(tprt<1e4&&tprt>0)
			{
				ep++;
				tprt*=10;
			}
			for(i0=1; i0<=ep;++i0)
			{
				strcat(ostr,"0");
				tprt/=10;
			}
			if(tprt==0)
			strcat(ostr,"0000");

		sprintf(var,"%i ",tprt);
		strcat(ostr,var);	

		for(i=2;i<=sf/6;i++)
		{
			ep=0;
			tprt=xlr.mant[i];
			while(tprt<1e5&&tprt>0)
			{
				ep++;
				tprt*=10;
			}
			for(i0=1; i0<=ep;++i0)
			{
				strcat(ostr,"0");
				tprt/=10;
			}
			if(tprt==0)
			strcat(ostr,"00000");

			sprintf(var,"%i ",tprt);
			strcat(ostr,var);
			if(i/10*10==i)
				strcat(ostr,"\n  ");
			if(i/100*100==i)
				{
				sprintf(var, " %i\n  ",(i*6)/100);
				strcat(ostr,var) ;  
			}

		}
		sprintf(var,"exp %i",xlr.exp);
		strcat(ostr,var);

} 

void invtan(NxlongR x,NxlongR* Itan,int dp)
{
	
	NxlongR xp,div,term,two;
	int sg=1,sdp=-dp;	
        
        if(x.exp<0)
            sdp=-dp+x.exp;
            
	

	*Itan=term=xp=x;
	itoxR(2,&two);
	itoxR(1,&div);
	
	prodxlr(x, x,&x,dp);
	while(term.exp>sdp)
	{
		sg*=-1;
		addXlr(div,two,&div,dp);
		prodxlr(xp, x,&xp,dp);
 		Divlxr( div,xp,&term,dp);
		term.sgn=sg;		//term.sgn *=sg
		addXlr(*Itan,term,Itan,dp);
	
	}
	return;
	
}


void CosSin(NxlongR x,int cs,int Ndp,NxlongR* CoS)	//cs 0,cos or 1,sin 
{
	NxlongR xn,fac,f2,t;
	int fct,sg;
	
	itoxR(1,CoS);
	if(cs==1)
	*CoS=x;
	xn=*CoS;
	
		prodxlr(x,x,&x,Ndp);
		fct=2;
		if(cs==1)
		fct=3;

		itoxR(2,&fac);
		if(cs==1)
		itoxR(6,&fac);

		f2.exp=0;
		sg=-1 * CoS->sgn;
		while (f2.exp>-Ndp)
		{
			
			prodxlr(x,xn,&xn,Ndp);
			
			Divlxr(fac,xn,&f2,Ndp);
			f2.sgn=sg;
			for(int i=1;i<=2;i++)
			{
				fct++;
				itoxR(fct,&t);
				prodxlr(fac,t,&fac,Ndp);
			}
			
			addXlr(f2,*CoS,CoS,Ndp);
			
			sg*=-1;
		}	
		return;	
}
void NRoot(NxlongR A, int N,NxlongR* x,int Ndp)  //05/22   /uses
// x(n+1)=ANx(n)/(AN-y) , y = x(n)^N -A
{
	NxlongR one,y,na;
	int i, i2=1,ye= A.exp,aN;	
	char istring[100];
	Ndp*=1.2;
	
	
	aN=N;
	itoxR(1,&one);
	if (N < 0)
	{
		printf("\nN is %d\n",N);
		Divlxr( A,one,&A,1.2*Ndp);
		crstrXlr(A,15,istring);
		printf("\n 1/A 	%s ",istring);
		aN *= -1;
	}
	
	itoxR(aN,&na);
	*x=y=one;
		             
    if((A.exp==-1 && A.mant[1]<800000) ||A.exp<-1 )
     {
		 Divlxr( na,A,&y,Ndp);
         if(A.exp<-5)
		{
		    y.exp=(A.exp-aN+1)/aN;
		} 
         
         *x=y;
	 }	
	prodxlr(*x, *x,&y,Ndp);
	for(i=3;i<=aN;++i)
	{
		prodxlr(*x, y,&y,Ndp);
	}	
	prodxlr(na, A,&na,Ndp); 
	A.sgn*=-1;       // -A
	
	
	addXlr(y,A,&y,Ndp);	
	
	addXlr(y,na,&y,Ndp);
	Divlxr( y,*x,&y,Ndp);	
		
	prodxlr(na, y,x,Ndp);	
	
	while(A.exp-ye< 0.9*Ndp )	
	{			
		prodxlr(*x, *x,&y,Ndp);		
		for(i=3;i<=aN;++i)
		{
			prodxlr(*x, y,&y,Ndp);
		}
		addXlr(y,A,&y,Ndp);	
				
		if (y.sgn == 0)
			break;
		ye = y.exp;
		
		addXlr(y,na,&y,Ndp);		
		
		Divlxr( y,*x,&y,Ndp);
				
		prodxlr(na, y,x,Ndp);		
		i2++;	             
	}
	if(N<0)
		printf("\nN is %d\n",N);
	return;
}			



// rdXlr requires global statements "namespace std;" and "ifstream InFile;"
/*void rdXlr(char fn[25],NxlongR* x,int mntit,int ofl)
{
	int i,mt;
	char istring[100];
	float m1;
	//InFile.open(fn, ios::in);
    for(i=0;i<ofl; i++)
    {
    	InFile.getline(istring,'\n');    	//InFile.getline(istring,100,'\n');
    }	
    
    	
	sscanf(istring,"%f %d %d %d %d %d %d %d %d %d",&m1, \
	&x->mant[2],&x->mant[3],&x->mant[4],&x->mant[5],&x->mant[6],&x->mant[7], \
	&x->mant[8],&x->mant[9],&x->mant[10]);
	
	x->sgn=1;
	if (m1<0)
		x->sgn=-1;
		
	mt=fabs(m1)*1e5;
	
	x->mant[1]=mt;
	
	
	for(i=1;i<mntit/10; i++)
	{	if (i/10*10==i)
			InFile.getline(istring,'\n'); 
			
		InFile.getline(istring,'\n');
		sscanf(istring,"%d %d %d %d %d %d %d %d %d %d",&x->mant[i*10+1], \
		&x->mant[i*10+2],&x->mant[i*10+3],&x->mant[i*10+4],&x->mant[i*10+5],&x->mant[i*10+6],&x->mant[i*10+7], \
		&x->mant[i*10+8],&x->mant[i*10+9],&x->mant[i*10+10]);
			
	}	
	InFile.getline(istring,'\n');
	InFile.getline(istring,'\n');
	sscanf(istring," exp %d",&x->exp);	
}
*/

