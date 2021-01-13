/*
		Here we define the structure describing our high precision 
		number and create the basic functions we require to manipulate and maintain them.
		Functions include the basic arithmetic functions of addition subtraction 
		multiplication and division ae well as some other useful process such as conversion 
		of an integer into a type of our number, direct calculation of a high precision 
		cosine or sine value or the high preecision calculation of a root.
*/


#include <math.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <time.h>
#include <cstring>

#define ncells 700		//remember you need extra cells for accuracy in division using Divlx

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

struct NxlongR
{
	int mant[ncells]; //each cell of array mant carries 6 digits of our extended real number
	int exp;	//exp is the exponent of our first non zero digit of mant[1]
	int sgn;	//overall sign -1 or 1 or if  NxlongR==0 sgn =0
};

//ifstream InFile;
FILE* InFile;

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

inline void xnegc(NxlongR* c,int dp) ///25/03/14
{
	double Dbl;
	int i,i2,ep,shft,shft2;
        for(i=ncells-1;i>(dp+5)/6+1;i--)
            c->mant[i]=0;
        
	for (i=(dp+5)/6;i>=2;i--)	// check for carry-up and neg results
		if(c->mant[i]<0)
		{
			c->mant[i-1]--;
			c->mant[i]+=1000000;
		}
		for (i=(dp+5)/6;i>1;i--)
		{
			if(c->mant[i]>=1e6)
			rippleup(c);	
		}
			
	if(c->mant[1]>=1e6)
	{
		Dbl=(double)c->mant[1];
		ep=0;

		while(Dbl>=1e6)
		{
			Dbl/=10;
			ep++;
		}
		c->mant[1]=Dbl;
		Dbl-=(double)c->mant[1];
		
		for(i=2;i<=(c->exp+dp)/6+1;i++)
		{
                    int j;
			for( j=1;j<=6+ep;j++)
				Dbl*=10;

			Dbl+=c->mant[i];
			for(j=1;j<=ep;j++)
				Dbl/=10;

			c->mant[i]=Dbl;
			Dbl-=(double)c->mant[i];
		}
		c->exp+=ep;
	}

	ep=0;
	i=1;
	if(c->mant[1]==0)
	{	
		i=2;
			while(c->mant[i]==0)
				i++;
		for(i2=1;i2<dp/6+2;i2++)
		{
			c->mant[0]=c->mant[i+i2-1];
			c->mant[i2]=c->mant[0];
		}

	
		if(i*6>=dp )
		{

			c->sgn=0;
			return;
		}
        if(c->mant[1]<0)
        {
            c->sgn*=-1;
            for(i2=1;i2<=(dp+5)/6;i2++)
            {
                c->mant[i2]*=-1;
            }
            
        }
	}
	
	while(c->mant[1]<1e5)
	{
		ep++;
		c->mant[1]*=10;
	}
	//if(i>1)
	c->exp-=ep+(i-1)*6;
	if(ep>0)
	{
		
			
		shft2=shft=c->mant[2];

		for( i2=1;i2<=6-ep;i2++)
			shft/=10;

		c->mant[1]+=shft;		

	
		for(i2=1;i2<=6-ep;i2++)
			shft*=10;

		shft=shft2-shft;


		i2=2;
		while((i2+5)*6<=dp)
		{
                    int i3;

			for( i3=1;i3<=ep;i3++)
				shft*=10;
			c->mant[i2]=shft;
			

			shft=shft2=c->mant[i+i2];

			for( i3=1;i3<=6-ep;i3++)
			shft/=10;
			c->mant[i2]+=shft;

			for( i3=1;i3<=6-ep;i3++)
			shft*=10;
			shft=shft2-shft;

			i2++;

		}
	}
		
}






// addXlr 19/03/14
void addXlr(NxlongR a, NxlongR b,NxlongR* c,int dp)
{
	NxlongR b2;
	int bign= a.exp-b.exp, xpb6=bign/6,i,i2,shft;
	//double Dbl;

	if(a.sgn==0)
	{
		*c=b;
		return;
	}
	if(b.sgn==0)
	{
		*c=a;
		return;
	}
        if(xpb6*xpb6>=(ncells-1)*(ncells-1))
        {
            if(bign>0)
            {
                *c=a;
                return;
            }
            *c=b;
            return;
        }

	for(i=1;i<ncells-1;i++)
		c->mant[i]=0;

	
	if(bign>0)
	{	
		for(i=1;i<=xpb6;i++)
		c->mant[i]=a.mant[i];

		int rem=bign-xpb6*6;
		b2.mant[1]=0;
		for(i=1;i<=dp/6+1;i++)		// b must be aligned before adding
		{
			shft=b.mant[i];
			for(i2=1;i2<=rem;i2++)
			shft/=10;

			b2.mant[i]+=shft;
			for(i2=1;i2<=rem;i2++)
			shft*=10;

			shft=b.mant[i]-shft;
			for(i2=1;i2<=6-rem;i2++)
			shft*=10;
			b2.mant[i+1]=shft;
				
		}
		for(i=dp/6+1;i>xpb6;i--)
						

		switch(a.sgn*b.sgn)
		{
			case	1:
			{
				c->mant[i]=a.mant[i]+b2.mant[i-xpb6];
				break;
			}
			
			default:
			{
				c->mant[i]=a.mant[i]-b2.mant[i-xpb6];
				break;
			}
	
		
		}
		c->exp=a.exp;
		c->sgn=a.sgn;
	}

	if(bign<0)
	{	
		bign=-bign;
		xpb6=bign/6;

		for(i=1;i<=xpb6;i++)
		c->mant[i]=b.mant[i];
		int rem=bign-xpb6*6;

		b2.mant[1]=0;
		for(i=1;i<=+dp/6+1;i++)		// a must be aligned before adding
		{

			shft=a.mant[i];
			for(i2=1;i2<=rem;i2++)
			shft/=10;

			b2.mant[i]+=shft;
			for(i2=1;i2<=rem;i2++)
			shft*=10;

			shft=a.mant[i]-shft;
			for(i2=1;i2<=6-rem;i2++)
			shft*=10;
			b2.mant[i+1]=shft;
			
	
		}
		for(i=dp/6+1;i>xpb6;i--)
						

		switch(a.sgn*b.sgn)
		{
			case 1:
			{
				c->mant[i]=b.mant[i]+b2.mant[i-xpb6];
				break;
			}

			default:
			{
				c->mant[i]=b.mant[i]-b2.mant[i-xpb6];
				break;
			}
		}

		c->exp=b.exp;
		c->sgn=b.sgn;
					
	}


	if(bign==0)
	{	
		i=0;
		if(a.sgn*b.sgn==1)
		while(i*6<=dp)
		{
			c->mant[i]= b.mant[i]+a.mant[i];
			++i;
		}
		c->sgn=b.sgn;
		c->exp=b.exp;

		if(a.sgn*b.sgn==-1)
		{	
			
			while(a.mant[i+1]==b.mant[i+1])
				i++;
			if(i*6>=dp)
			{
				c->mant[1]=0;
				c->sgn=0;
				return;
			}

			i2=1;
			while(i+i2<=ncells-1)
			{
			
					{
						if(a.mant[i+1]>b.mant[i+1])					
						c->mant[i2]= a.mant[i+i2]-b.mant[i+i2];

						else				
						c->mant[i2]= b.mant[i+i2]-a.mant[i+i2];					
					
					}
								
				i2++;
			}
			c->exp=a.exp-i*6;

		
				if(b.mant[i+1]>a.mant[i+1])
				c->sgn=b.sgn;
				else
				c->sgn=a.sgn;
			
		}
		if(b.mant[1]>a.mant[1])
		{	
			 i2=1;
			while(i2*6<=dp)
			{
				switch(a.sgn*b.sgn)
				{
					case 1:
					{
						c->mant[i2]= a.mant[i2]+b.mant[i2];
						break;
					}

					default:
					{
						c->mant[i2]= b.mant[i2]-a.mant[i2];
						break;
					}
				}
				i2++;
			}
			c->exp=b.exp;
			c->sgn=b.sgn;
		}



		if(b.mant[1]<a.mant[1])
		{	
			 i2=1;
			while(i2*6<dp)
			{
			switch(a.sgn*b.sgn)
				{
					case 1:
					{
						c->mant[i2]= a.mant[i2]+b.mant[i2];
						break;
					}

					default:
					{
						c->mant[i2]= a.mant[i2]-b.mant[i2];
						break;
					}
				}	
				i2++;
			}
			c->exp=a.exp;
			c->sgn=a.sgn;
		}
		
	}
	
	xnegc(c,dp);
   xnegc(c,dp);
			
	if (c->mant[1]==0)
		c->sgn=0;

}

inline void prodxlr(NxlongR a, NxlongR b,NxlongR* c,int dp) //25/03/14
{
	int i,j,sf1b6,sf2b6,f,m,bk,ep=0;	//sufx=0;

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

//09/01/14
void Divlxr(NxlongR divisor,NxlongR divend,NxlongR* quot,int dp)
{
	
	NxlongR rem,Int2,quot2;	//Int1,
	//char str[100];
	double Dvisor,Dvend,Dquot;
	int ex;
	dp*=1.3;

	if(divisor.sgn==0)
	{
		quot->exp=12345;
		quot->sgn=divisor.sgn*divend.sgn;
		return;
	}
	if(divend.sgn==0)
	{
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

	
	int i=1,sf=12;			//sfc, i2,   sf for quotient is different from remainder may require increasing for validity

	
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

	prodxlr(*quot,divisor,&Int2,dp);

	divend.sgn=1;
	Int2.sgn=-1;
		
	addXlr(divend,Int2,&rem,dp*0.75);

	if(rem.sgn==0)
		return;

		
	while(divend.exp-rem.exp<0.8*dp-1)		//dp =1.3 dp(in)
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

		
		/*	//diagnostic block
			if(quot2.exp<-40)
			{
				crstrXlr(divend,8,str);
				fprintf(OutFile,"\n\n\r %s\n\r",str);
				crstrXlr(quot2,15,str);
				fprintf(OutFile,"\n\r %s\n\r",str);
				crstrXlr(rem,15,str);
				fprintf(OutFile,"\n\r %s\n\r",str);
			}
//				*/
	
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
			while(i<ncells-1)
			{	
				i++;		
				quot2.mant[i]=0;
				
			}
			xnegc(&quot2,dp);

						
				addXlr(*quot,quot2,quot,dp*0.8);

				prodxlr(*quot,divisor,&Int2,dp);
				Int2.sgn*=-1;
				addXlr(divend,Int2,&rem,dp*0.75);
	
		

			if(rem.sgn==0)
				return;		
		
	}   // end while
		
//	prodxlr(*quot,divisor,&Int1,1.25*dp);		//diagnostic

}
void itoxR(int intg,NxlongR* xlo)
{
	double tmp1;
	int ep=0,i;//,Ndp=(ncells+1)/6
	
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
	tprt= t2=xlr.sgn*xlr.mant[1];
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
			if (i / 100* 100== i)
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
		sg=-1;
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

void interp(NxlongR h,NxlongR* x1,NxlongR A,int N,int dp)
{
	NxlongR f1,f2,f3,x2,x3;	//,hlf;
	int m1h[11]={0,110,109,108,107,106,105,104,103,102,101},i;	//stop cycleing with h same

	dp*=1.5;
	/*hlf.mant[1]=500000;
	for(i=3;i<=ncells-1;++i)
	hlf.mant[i]=0;

	hlf.exp=-1;
	hlf.sgn=1;*/

	
	//A.sgn*=-1;   //A is already -ve
	

		prodxlr(*x1, *x1,&f1,dp);

		for(i=3;i<=N;++i)
		{
	
			prodxlr(*x1, f1,&f1,dp);
			
		}
		addXlr(f1,A,&f1,dp);
	

		addXlr(*x1,h,&x2,dp);
		prodxlr(x2, x2,&f2,dp);


		for(i=3;i<=N;++i)
		{
			prodxlr(x2, f2,&f2,dp);
			
		}
		addXlr(f2,A,&f2,dp);

		f1.sgn*=-1;
		addXlr(f1,f2,&f3,dp);		//	f2-f1
		
		if (f3.sgn==0)
		{
			return;
			
		}
		else
		Divlxr( f3,f1,&x3,dp);

		f1.sgn*=-1;


		prodxlr(x3,h,&h,1.2*dp);
	
	while(h.exp-(A.exp/N)>-0.68*dp)
	//while(h.exp-(A.exp/N)>-dp)
	{
		for(i=1;i<=9;++i)
		m1h[i]=m1h[i+1];

		m1h[10]=h.exp;

		if(m1h[10]>=m1h[1])
		{
		//	itoxR(-1,&hlf);
		//	prodxlr(h,hlf,&h,dp);
			h.sgn*=-1;
			addXlr(*x1,h,x1,dp);
			return;
		}


		addXlr(*x1,h,&x3,dp);
		
		
	
		prodxlr(x3, x3,&f3,dp);
		for(i=3;i<=N;++i)
		{
			prodxlr(x3, f3,&f3,dp);
			
		}
		
		addXlr(f3,A,&f3,dp);
		if(f3.sgn==0)
		{
			*x1=x3;
			return;
		}
		if(f2.exp<f1.exp)
		{
			*x1=x2;
			x2.sgn*=-1;
			addXlr(x2,x3,&h,dp);

		}

			prodxlr(*x1, *x1,&f1,dp);

		for(i=3;i<=N;++i)
		{
	
			prodxlr(*x1, f1,&f1,dp);
			
		}
		addXlr(f1,A,&f1,dp);
	

		addXlr(*x1,h,&x2,dp);
		prodxlr(x2, x2,&f2,dp);


		for(i=3;i<=N;++i)
		{
			prodxlr(x2, f2,&f2,dp);
			
		}
		addXlr(f2,A,&f2,dp);

		f1.sgn*=-1;
		addXlr(f1,f2,&f3,dp);		//	f2-f1
		
		if (f3.sgn==0)
		{
			return;
			
		}
		else
		Divlxr( f3,f1,&x3,dp);

		f1.sgn*=-1;


		prodxlr(x3,h,&h,1.2*dp);

		
	}

	


}


void NRootx(NxlongR A, NxlongR h,int N,NxlongR* x,int Ndp) // revised 09/20
{
	NxlongR f,fda,n;
	int i,c1=1,hxp; //ep=0,expo=0,sgh=2,sg=1,double xdbl;

			
	/*/double Adbl= A.mant[1]/1e5*pow(10,A.exp);
	itoxR(1,x);
	for(i=4;i<=ncells-1;i++)
	x->mant[i]=0;
	
	
	h.exp=x->exp=(A.exp+N-1)/N;*/
	hxp=h.exp;
		
	
	A.sgn=-1;
	while (h.exp-x->exp>-Ndp)//(h.exp-(A.exp/N)>-Ndp)
	{
		if(c1>2)
		if(hxp<h.exp)
		{
	
			
				h.sgn*=-1;

				addXlr(*x,h,x,Ndp);
				
				interp(h,x,A,N,Ndp);
				
				prodxlr(*x, *x,&f,Ndp);

				for(i=3;i<=N;++i)
					prodxlr(*x, f,&f,Ndp);
					
				addXlr(f,A,&f,Ndp);
				return;


		}
		c1++;
	
		prodxlr(*x, *x,&f,Ndp);

		for(i=3;i<=N;++i)
		{
			if(i==N)			//counts as N-1 
			fda=f;

			prodxlr(*x, f,&f,Ndp);
			
		}

		if (N==2)
		fda=*x;

		//f.sgn*=1;
		
		addXlr(f,A,&f,Ndp);
		//A.sgn*=-1;
		if(f.sgn==0)
			return;
		
       
		itoxR(N,&n);
		//prodxlr(n, fda,&fda,1.4*Ndp);
		prodxlr(n, fda,&fda,1.2*Ndp);
		
		
		hxp=h.exp;
		cout <<"\r  "<<hxp;

		//Divlxr( fda,f,&h,Ndp);
		Divlxr( fda,f,&h,200);
		for(i=200/6+1;i<=Ndp/6+1;++i)
			h.mant[i]=0;

				//sgh=f.sgn;
			
				h.sgn*=-1;

				addXlr(*x,h,x,Ndp);
	}


}
