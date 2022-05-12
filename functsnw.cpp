
/*
		Here we define the structure describing our high precision 
		number and create the basic functions we require to manipulate and maintain them.
		Functions include the basic arithmetic functions of addition subtraction 
		multiplication and division ae well as some other useful process such as conversion 
		of an integer into a type of our number, direct calculation of a high precision 
		cosine or sine value or the high preecision calculation of a root.
*/

/*
A compilation of the most updated functions of ArithXlr
includes corrected addXlr() , Divlxr(), NRoot() with its adjunt interp()
Basic functions  CosSin(), invtan() and itoxR() are included
also the print function crtxlxr() in this form it displayed dp/100 in the line between blocks


inline void addXlr(NxlongR a, NxlongR b, NxlongR * c, int dp)
inline void prodxlr(NxlongR a, NxlongR b,NxlongR* c,int dp) 
inline void Divlxr(NxlongR divisor, NxlongR divend, NxlongR * quot, int dp)
void itoxR(int intg,NxlongR* xlo)
void invtan(NxlongR x,NxlongR* Itan,int dp)
inline void crstrXlr(NxlongR xlr,int sf,char* ostr)
void CosSin(NxlongR x,int cs,int Ndp,NxlongR* CoS)	//cs 0,cos or 1,sin 
void NRoot(NxlongR A, int N,NxlongR* x,int Ndp) // revised improved 28/03/21  

*/

	
	# include <stdlib.h>
	# include <math.h>
	# include <stdio.h>
	# include <string.h>
	# include <iostream>	
	# include <time.h>
	# define ncells 2000 //remember you need extra cells for accuracy in division using Divlx

/*Developed by Olatunde Adeyemo �2013
PIsf3 is derived from expo.cpp to calculate exp(x)

let us assume we are working to n dp through int Ndp
*/
//#define USEC_TO_SEC 1.0e-6

using namespace std;



struct NxlongR {
	int mant[ncells]; //each cell of array mant carries 6 digits of our extended real number
	int exp; //exp is the exponent of our first non zero digit of mant[1]
	int sgn; //overall sign
};

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
	//if(i>1)
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



			
			
//True working ver of Divlxr as at 21/12/13
//09/01/14  //corrected 25/03/21 /03/05/21
inline void Divlxr(NxlongR divisor, NxlongR divend, NxlongR * quot, int dp) {

	NxlongR rem,  Int2, quot2,one, X;//Int1,
        dp = 1.4*dp;
	double Dvisor, Dvend, Dquot;
	int ex,i;
	
	i=1;
	one.mant[i]=100000;
	while (i < ncells - 1) {
		i++;
		one.mant[i] = 0;
	}
	one.exp= 0;
	one.sgn= 1;
	

	if (divisor.sgn == 0) {
		quot -> exp = 12345;
		quot -> sgn = divisor.sgn * divend.sgn;
		return;
	}
	if (divend.sgn == 0) {
		quot -> sgn = 0;
		return;
	}

	Dvisor = divisor.mant[1] * 1e-5 + divisor.mant[2] * 1e-11 + divisor.mant[3] * 1e-17 ;
	Dvend = divend.mant[1] * 1e-5 + divend.mant[2] * 1e-11 + divend.mant[3] * 1e-17 ;
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
	quot -> exp = ex + divend.exp - divisor.exp;
	quot -> sgn = divend.sgn * divisor.sgn;
	
	
	Dquot *= 1e5;
	quot -> mant[1] = Dquot;


	int  sf = 12; //, sfc, i2, m10sf for quotient is different from remainder may require increasing for validity
	i = 1;

	while ((i + 1) <= (sf + 5) / 6) {
		Dquot -= quot -> mant[i];
		i++;
		Dquot *= 1e6;
		quot -> mant[i] = Dquot;
	}
	while (i < ncells - 1) {
		i++;
		quot -> mant[i] = 0;

	}

	xnegc(quot, dp);
        
        Int2 = one;
        Int2.sgn = -quot->sgn;
        addXlr(Int2, *quot, &Int2, dp);
        if(Int2.exp - divend.exp <= -12)  
        {            
           
            divisor.sgn *= -1;
            addXlr(divend, divisor, & X, dp);
            divisor.sgn *= -1;

            Dvend = X.mant[1] * 1e-5 + X.mant[2] * 1e-11 + X.mant[3] * 1e-17;
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
            rem.exp = ex + X.exp - divisor.exp;
            rem.sgn = divisor.sgn * X.sgn; 

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
        }
               
        prodxlr( * quot, divisor, & Int2, dp);
	
	Int2.sgn *= -1;		//******

	addXlr(divend, Int2, & rem, dp);

	if (rem.sgn == 0)
		return;

	while (divend.exp - rem.exp <=  dp*0.8 ) //dp =1.5 dp(in)0.68,
	{   
            Dvend = rem.mant[1] * 1e-5 + rem.mant[2] * 1e-11+ rem.mant[3] * 1e-17;
            Dvisor = divisor.mant[1] * 1e-5 + divisor.mant[2] * 1e-11+ divisor.mant[3] * 1e-17;
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

                quot2.exp = ex + rem.exp - divisor.exp;
                quot2.sgn = rem.sgn * divisor.sgn;


                Dquot *= 1e5;
                quot2.mant[1] = Dquot;


                i = 1;
                while ((i + 1) <= (sf + 5) / 6) {
                        Dquot -= quot2.mant[i];
                        i++;
                        Dquot *= 1e6;
                        quot2.mant[i] = Dquot;
                }
                while (i < ncells - 1) {
                        i++;
                        quot2.mant[i] = 0;

                }                

                addXlr( * quot, quot2, quot, dp);
         
            quot -> sgn = divend.sgn * divisor.sgn;
            prodxlr( * quot, divisor, & Int2, dp);
            Int2.sgn *= -1;
           
            addXlr(divend, Int2, & rem, dp);           
            if (rem.sgn == 0)
                    return;
	}       
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
		



inline void crstrXlr(NxlongR xlr,int sf,char* ostr)
{
	char var[30];
	int i,tprt,t2,i0,ep=0;
	tprt= t2=xlr.sgn*xlr.mant[1],i=1;
	for(i=1;i<=5;i++)
		tprt/=10;

	sprintf(ostr,"\r  %i.",tprt);
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
			//addXlr(t,one,&t,Ndp);
		
			prodxlr(x,xn,&xn,Ndp);
			
			Divlxr(fac,xn,&f2,Ndp);
			f2.sgn=sg;
			for(int i=1;i<=2;i++)
			{
				fct++;
				itoxR(fct,&t);
				prodxlr(fac,t,&fac,Ndp);
			}
			
		//	if(cs==0)
		//	prodxlr(x,xn,&xn,Ndp);

			addXlr(f2,*CoS,CoS,Ndp);
			sg*=-1;
		}
		//cos,sin =f;
		return;	
}


void NRoot(NxlongR A, int N,NxlongR* x,int Ndp) // revised 28/03/21
{												// using a(i+1) = a(i)/(1+(a(i)^n - A)/nA)	
	NxlongR f,h,fda,n,mone,xo;					// as in fnsconvrt.cpp
	int i;
	char istring[100];


	fda=A;
	fda.sgn *=-1;
	
	h.exp =(A.exp+N-1)/N;	
	itoxR(1,&mone);		
	*x =mone;
	mone.sgn=-1;	
	itoxR(N,&n);	
	//A.sgn=-1;
	while (h.exp-(A.exp/N)>-Ndp) // h has been generated	
	{
		//printf("\n\rh.exp-(A.exp/N) %d\n",h.exp-(A.exp/N));
		
		xo=*x;
		xo.sgn *=-1;
	
		prodxlr(*x, *x,&f,1.2*Ndp);

		for(i=3;i<=N;++i)
		{
			
			prodxlr(*x, f,&f,1.2*Ndp);
			
		}	
		addXlr(f,fda,&h,1.2*Ndp);
		//crstrXlr(h,15,istring);
		//printf("\r\nrem %s\n\r",istring);

		if(h.sgn==0)
			return;
			
		Divlxr( A,f,&f,1.2*Ndp);
		
		addXlr(f,mone,&f,1.2*Ndp);
		addXlr(f,n,&f,1.2*Ndp);
		
		Divlxr( f,n,&f,1.2*Ndp);	
		
		prodxlr(f, *x,x,1.2*Ndp);
		addXlr(*x,xo,&h,Ndp*1.2);
		//crstrXlr(h,15,istring);
		//printf("\r\nh %s\n\r",istring);
		//fprintf(OutFile,"\r\nh %s\n\r",istring);		
	}	

}
