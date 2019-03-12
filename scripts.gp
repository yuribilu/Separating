/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



PARI scripts used for the article 

"Separating singular moduli and the primitive element problem"

Yu. Bilu, B. Faye, H. Zhu

Table of contents

A. General scripts

minabs(v) smallest  absolute value in a vector
distance_in(v) minimal distance between entries of a vector
distance(v,w) minimal distance between entries of two vectors
areprop(f,g) are two polynomials proportional

B. Scripts on j,  dicriminants and singular moduli

jprime(z) = j'(z)
jsecond(z)=j''(z)
taus(Delta) finding the tau_x of given discriminant Delta
sigmods(Delta) finding the singular moduli of given discriminant Delta
deltas(X) discriminants up to X
deltas_ge_h(X,h) list of discriminants up to X with class number at least h
mindistconj(X) smallest distance between two distinct conjugate singular moduli up to disriminant X
sigmodsupto(X) list of the singular moduli up to disriminant X
mindist(X) smallest distance between two distinct singular moduli up to disriminant X
separeal(v), sepanonreal(v) separating (non-)totally real discriminants

C. Scripts for the article

cor53j(X), cor53j1728(X), cor53jprime(X) verifications for Corollary 5.3
theorem61(X) verification for Theorem 6.1
baddeltas(N) discriminants satisfying (8.4) or (8.5)
simag(v) smallest absolute value of imaginary parts of (v[1]-v[i])/(v[j]-v[k])
smallest_simag(v) smallest simag for a list of discriminants
aredifsprop(v) proportional differences of polynomials
iscrossratiorational(f) no rational cross-ratio of conjugates of an algebraic number generating a Galois extension of Q
ratrootdifquot(f,g) can difference of roots of of f divided by that of g be rational?
sec83(v) verification for Section 8.3 
list83 list for Section 8.3




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/ 

/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

A. General scripts

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/


/*--------------------------------------------------------------------------------
smallest  absolute value in a vector
input: vector with real or complex entries
output: smallest absolute value of its entries
-----------------------------------------------------------------------------------*/

minabs(v)=
{
my(w,n);
n=length(v);
w=vector(n,i, abs(v[i]));
return(vecmin(w));
}

/*-------------------------------------------------------------------------------
minimal distance between entries of a vector
input: vector v with real or complex entries
output: min|v[i]-v[j]|, i\ne j
------------------------------------------------------------------------------------*/

distance_in(v) = 
{ 
my(i,j,n,d);

n=length(v);

if (n==1, print ("Vector must have at least 2 entries"); 
return(0););

d=abs(v[1]-v[2]);

for (i=1,n-1,
for (j=i+1,n,
d=min(d,abs(v[i]-v[j]));););
return(d);
}

/*-------------------------------------------------------------------------------
minimal distance between entries of two vectors
input: vectors v,w with real or complex entries
output: min|v[i]-w[j]| 
------------------------------------------------------------------------------------*/

distance(v,w) = 
{ 
my(m,n,d);

m=length(v);
n=length(w);



/*d=abs(v[1]-w[1]);*/

d=+oo;

for (i=1,m,
for (j=1,n,
d=min(d,abs(v[i]-w[j]));););
return(d);
} 

/*-------------------------------
are two polynomials proportional
input: non-zero polynomials f,g in the same main variable
output: 1 if f/g is constant; 0 if not
----------------------------------*/

areprop(f,g)=
{
my(h);
if(f*g==0, print("f or g is 0"); return(0););
if(poldegree(f)!=poldegree(g),return(0););
h=pollead(f)*g-pollead(g)*f;
if(h==0, return(1););
return(0); 
}

/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

B. Scripts on j, discriminants and singular moduli

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/


/*-----------------------------------------------
j'(z)
--------------------------------------------------*/

jprime(z) =
{ my(j,L,E4,E6);
  L = [z,1]; 
  E4 = elleisnum(L,4);
  E6 = elleisnum(L,6);
  j = 1728*E4^3/(E4^3 - E6^2);
  - E6/ (E4*2*I*Pi) * j;
}


/*-----------------------------------------------------------
j''(z)
---------------------------------------------------------*/

jsecond(z) =
{ my(j,L,E4,E6);
  L = [z,1]; 
  E2 = elleisnum(L,2);
  E4 = elleisnum(L,4);
  E6 = elleisnum(L,6);
  2^5*3^2*(3*E4^4+4*E4*E6^2-E2*E4^2*E6)/((E4^3 - E6^2)*(2*Pi*I)^2);
}



/*--------------------------------------------------------------------------
finding the tau_x of given discriminant Delta
input: Delta (must be <0, =0,1mod4, otherwise error message)
output list of all \tau_x
----------------------------------------------------------------------------- */

taus(Delta) =
{
    my(a,b,c,clno,m,v);

     if (Delta>0, 
    print("Bad Delta");
    return(0););


    if (Delta%4==2, 
    print("Bad Delta");
    return(0););
	
	    if (Delta%4==3, 
    print("Bad Delta");
    return(0););
    
   
    clno  =qfbclassno(Delta);

    X = abs(Delta);
	
	v=vector(clno,i,0); 
	
    maxa = floor(sqrt(X/3));



    for(a=1,maxa,
	for(b=-a+1,a,
	    if((b^2-Delta)%(4*a) != 0, next;);
	    c = (b^2-Delta)/(4*a);
	    
	    /*
	      Note that b>-a, next two ifs check that
	      we are in the fundamental domain.
	    */	 

	    if(a>c,next;);
	    if( (a==c) && (b<0),next;);

	    /* (a,b,c) must be coprime. */
	    if(gcd(a,gcd(b,c))!=1,next;);
	    
	    
	    m = m + 1;
		
			
v[m]=(b+I*sqrt(X))/(2*a);		 
	);
    );

/* Sanity check. Make sure our n matches with Pari's computation */
/* Used only for debugging purposes */
    if(m != clno,
        return("FAIL (internal error)");
    	return(0);
    ); 
    return(v);
}


/*------------------------------------------------------------
finding the singular moduli of given discriminant Delta
--------------------------------------------------------------*/

sigmods(Delta) = 
{
    my(v,w,clno);

     if (Delta>0, 
    print("Bad Delta");
    return(0););


    if (Delta%4==2, 
    print("Bad Delta");
    return(0););
	
	    if (Delta%4==3, 
    print("Bad Delta");
    return(0););
    
   
    clno  =qfbclassno(Delta);
	
	w=taus(Delta);
	
	v=vector(clno,i,ellj(w[i])); 
	
    
    return(v);
}



/*--------------------------------------------------------------------
list of discriminants up to X
input: X\ge 3
output: discriminants Delta  in [-X,0]
----------------------------------------------------------------------*/

deltas(X)=
{
my(m,n,i);
if (X<3, print("X too small");return(0););
m=floor((X+1)/4);
n=floor(X/4);
v=vector(m+n,i,0);
for(i=1,m,
v[2*i-1]=1-4*i;);
for(i=1,n,
v[2*i]=-4*i;);
return(v);
}

/*--------------------------------------------------------------------
list of discriminants up to X with class number at least h
input: X\ge 3, h\ge 1
output: discriminants Delta  in [-X,0] with h(Delta) \ge h
----------------------------------------------------------------------*/

deltas_ge_h(X,h)=
{
my(u,v,w,n,i,j);
if (X<3, print("X too small");return(0););
if (h<1, print("h too small");return(0););
u=deltas(X);
n=length(u);
v=vector(n,i,0);
j=1;
for(i=1,n,
if(qfbclassno(u[i])>=h, 
v[j]=u[i];
j=j+1;););
if(j=1, print("X too small");return(0););
w=vector(j-1,i,0);
w=v[1..j-1];
return(w);
}



/*------------------------------------------------------------------
smallest distance between distinct conjugate singular moduli up to disriminant X
input: X\ge 15
output: smallest distance 
Remark: used to prove Proposition 6.2
-----------------------------------------------------------------------*/

mindistconj(X)=
{
my(n,Delta,d,v);

if (X<15, print("X too small"); return(0););

d=distance_in(sigmods(-15));

for (n=15,floor(X), 
if (n%4==1,next;);
if (n%4==2,next;);
Delta=-n;
if (qfbclassno(Delta)==1,next;);
v=sigmods(Delta);
d=min(d,distance_in(v)););
return(d);
}


/*------------------------------------------------------------------
list of singular moduli of  disriminants bounded by X
input: X\ge 3 
output: list of singular moduli of discriminant bounded by X
-----------------------------------------------------------------------*/

sigmodsupto(X)=
{
my(a,b,c,maxa,maxc,maxcount,m,v,w,d);

 if (X<3, print("X too small"); return(0););

maxa=floor((X/3)^1/2);

maxcount=floor(maxa*(X+1)/2); 

/*
We have maxa possibilities for a;
for every given a we have 2a possibilities for b;
for every given (a,b) we have at most 
(X+b^2)/(4a) - a +1 possibilities for c.
We have clearly 
(X+b^2)/(4a) - a +1 \le (X+1)/(4a). 
Hence we have at most maxa*(X+1)/2 possible (a,b,c). 
*/
	
	v=vector(maxcount,i,0); 
	

m=0;


    for(a=1,maxa,
	for(b=-a+1,a,
maxc=floor((X+b^2)/(4*a));
if(maxc<a, next;);

for(c=a,maxc,	

	    
	    /*
	      Note that b>-a, next two ifs check that
	      we are in the fundamental domain.
	    */	 

	    if( (a==c) && (b<0),next;);

	    /* (a,b,c) must be coprime. */
	    if(gcd(a,gcd(b,c))!=1,next;);

Delta=b^2-4*a*c;

/* print(Delta); */

m=m+1;

	    
v[m]=ellj((b+(-Delta)^(1/2)*I)/(2*a));	  

/* print (v[m]); */

		 
	);
    );
	);
	
/*	print(m); */
	

	
	w=vector(m,i,0); 

  w=v[1..m];
 

 
    return(w); 
}



/*------------------------------------------------------------------
smallest distance between two distinct singular moduli x,y of  disriminants bounded by X
input: X\ge 4 
output: smallest distance 
Remark: used to prove Proposition 6.2
-----------------------------------------------------------------------*/

mindist(X)=
{
my(v,d);

 if (X<4, print("X too small"); return(0););

v=sigmodsupto(X);
 
 d=distance_in(v); 
 
    return(d); 
}

/*----------------------------------------------
separating totally real discriminants
input: v list of discriminants
output
sublist of v containing the discriminants with only real singular moduli
------------------------------------------------*/
separeal(v)=
{
my(i,j,n,u,w,f,Delta,h);
n=length(v);
w=vector(n,i,0);
j=1;
for(i=1,n,
Delta=v[i];
f=polclass(Delta);
h=poldegree(f);
r=length(polrootsreal(f));
if(h==r, 
w[j]=Delta;
j=j+1;););
u=w[1..j-1];
return(u);
}

/*-----------------------------------------
separating non-totally real discriminants
input: v list of discriminants
output
sublist of v containing the discriminants with at least one non-real singular modulus
--------------------------------------------*/
sepanonreal(v)=
{
my(i,j,n,u,w,f,Delta,h);
n=length(v);
w=vector(n,i,0);
j=1;
for(i=1,n,
Delta=v[i];
f=polclass(Delta);
h=poldegree(f);
r=length(polrootsreal(f));
if(h!=r, 
w[j]=Delta;
j=j+1;););
u=w[1..j-1];
return(u);
}

/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

C. Scripts for the article

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/



/*----------------------------------------------------------
verifications for Corollary 5.3
input: X\ge 7
output: 1 if (6.4),respectively (6.5), respectively (6.6) holds for all imaginary quadratic tau with discriminants \le X; 0 if not (together with the list of exceptions) 
-----------------------------------------------------------*/

cor53j(X)=
{
my(n,u,v,w,clno,Delta,d,rhs,tf);
if(X<7, print ("X too small");return(0));
tf=1;
u=deltas(X);
n=length(u);
for(j=3,n,
Delta=u[j];
clno=qfbclassno(Delta);
v=taus(Delta);
w=vector(clno,i,abs(ellj(v[i])));
rhs=700*Delta^(-3);
d=vecmin(w);
if(d<rhs,
print("bad discriminant " Delta);
tf=0;););
return(tf);
}



cor53j1728(X)=
{
my(n,u,v,w,clno,Delta,d,rhs,tf);
if(X<7, print ("X too small");return(0));
tf=1;
u=deltas(X);
n=length(u);
for(j=3,n,
Delta=u[j];
clno=qfbclassno(Delta);
v=taus(Delta);
w=vector(clno,i,abs(ellj(v[i])-1728));
rhs=2000*Delta^(-2);
d=vecmin(w);
if(d<rhs,
print("bad discriminant " Delta);
tf=0;););
return(tf);
}


cor53jprime(X)=
{
my(n,u,v,w,clno,Delta,d,rhs,tf);
if(X<7, print ("X too small");return(0));
tf=1;
u=deltas(X);
n=length(u);
for(j=3,n,
Delta=u[j];
clno=qfbclassno(Delta);
v=taus(Delta);
w=vector(clno,i,abs(jprime(v[i])));
rhs=40000*Delta^(-2);
d=vecmin(w);
if(d<rhs,
print("bad discriminant " Delta);
tf=0;););
return(tf);
}


/*---------------------------------------------
right-hand side of (6.1)
input u\ge v>0
-------------------------------------------*/

rhs61(u,v)=vecmin([800*v^(-4),20000*u^(-1)*v^(-3),700*u^(-3)])





/*-------------------------------------------------------------
verification for Theorem 6.1
input: X\ge 4, Y\ge 3
output: 1 if (6,1) holds for all singular moduli x, y with discriminants \le X, \le Y respectively; 0 if not (together with the list of exceptions) 
--------------------------------------------------------------------*/

theorem61(X,Y)=
{
my(u,v,m,n,d,rhs,tf,clno);
if(X<4, print ("X too small");return(0));
if(Y<3, print ("Y too small");return(0));
u=deltas(X);
v=deltas(Y);
m=length(u);
n=length(v);
tf=1;
for(i=2,m,
for(j=1,min(n,i-1),
rhs=rhs61(u[i],v[j]);
d=distance(sigmods(u[i]),sigmods(v[j]));
if(d<rhs,
print("bad discriminants " u[i] "," v[j]);
tf=0;);););
for(i=1,min(m,n),
rhs=rhs61(u[i],u[i]);
clno=qfbclassno(u[i]);
if(clno==1,next;);
d=distance_in(sigmods(u[i]));
if(d<rhs,
print("bad discriminants " u[i] "," u[i]);
tf=0;););
return(tf);
}


/*-----------------------------------------------
Discriminants satisfying (8.4) or (8.5)
input N\ge 20
output
all discriminants satisfying (8.4) or (8.5) not exceeding 2*N in absolute value
Remark: since we know that h(Delta)>6 when |Delta| > 3763, it is enough to run with N=1882
----------------------------------------------------*/

baddeltas(N)=
{
my(n,i,v,w,Delta,clno);
if (N<20, print("N too small"); return(0););
v=vector(N,i,0);
i=1;
for (n=39,2*N, 
if (n%4==1,next;);
if (n%4==2,next;);
Delta=-n;
clno=qfbclassno(Delta);
if ((clno==4||clno==5||clno==6)&&Delta%8==1, 
v[i]=Delta;
i=i+1;
next;);
if (clno==4&&(Delta%16==8||Delta%16==12),
v[i]=Delta;
i=i+1;
next;);
);
print(i-1 " bad Deltas");
return(v[1..i-1]);
}

/*------------------------------------------------------
smallest absolute value of imaginary parts of (v[1]-v[i])/(v[j]-v[k])
input: v vector of n distinct complex numbers
output: smallest imaginary parts of (v[1]-v[i])/(v[j]-v[k]) with 
1< i,j <n , j<k\le n
--------------------------------------------------------*/ 

simag(v)=
{ 
my(n,ima);

n=length(v);

if (n<3, print ("Vector must have at least 3 coordinates"); 
return(0););

ima=abs(imag((v[1]-v[2])/(v[2]-v[3])));

for (i=2,n-1,
for (j=2,n-1,
for (k=j+1,n,
ima=min(ima,abs(imag((v[1]-v[i])/(v[j]-v[k])))););););
return(ima);
}


/*------------------------------------
smallest simag for a list of discriminants
input: v a list of discriminants
output: computing simag for the singular moduli of each of the discrimints, and finding the smallest of them
---------------------------------------*/

smallest_simag(v)=
{
my(w,i,n,ima);
n=length(v);
w=vector(n,i, simag(sigmods(v[i])));
print("smallest simag in the list="vecmin(w));
return(w);
}



/*---------------------------------------
proportional differences of polynomials
input: vector of  n\ge 3 polynomials
output: 0 if for some  i,j,k with 1< i,j <n , j<k\le n the polynomials 
v[1]-v[i] and v[j]-v[k] are proportional; 1 if not
used only for debugging iscrossratiorational
----------------------------------------*/



aredifsprop(v)=
{
my(n);

n=length(v);

if (n<3, print ("there must be at least 3 polynomials"); 
return(0););
for (i=2,n-1,
for (j=2,n-1,
for (k=j+1,n,
if(areprop(v[1]-v[i],v[j]-v[k])==1, print("("v[1]")-("v[i]") and ("v[j]")-("v[k]") are proportional");return(0););
);););
return(1);
}


/*--------------------------------------
no rational cross-ratio of conjugates of an algebraic number generating a Galois extension of Q
input: a Q-irreducible monic polynomial f in Z[x] of degree at least 3 whose root generates a Galois extension of Q;
output: 0 if some cross-ratio of the roots is (non-trivially) rational; 1 if none.  
--------------------------------------*/



iscrossratiorational(f)=
{
my(v,n,K);
n=poldegree(f);
if (n<3, print("deg f<3"); return(0););
if(polisirreducible(f)==0, print("f is reducible"); return(0););
K=nfinit(f);
v=nfgaloisconj(K);
if(length(v)<n,
print(f "not Galois"); 
return(0););

for (i=2,n-1,
for (j=2,n-1,
for (k=j+1,n,

if(areprop(v[1]-v[i],v[j]-v[k])==1, print(v[1]","v[i]","v[j]","v[k]" have rational cross-ratio");return(0););

);););
return(1);
}



/*------------------------------------
no rational cross-ratio in a list of 2-elementary discriminants
input: list u of of 2-elementary discriminants with class number at least 4
output: 0 if some cross-ratio of the the singular moduli of some discriminant in the list  is (non-trivially) rational; 1 if none. 
--------------------------------------*/



norationalcrossratio(u)=
{
my(v,m,Delta,n,K);
m=length(u);
for(ell=1,m,
Delta=u[ell];
f=polclass(Delta);
n=poldegree(f);
if (n<4, print(Delta" has class number <4"); return(0););
K=nfinit(f);
v=nfgaloisconj(K);
if(length(v)<n,
print(Delta "not 2-elementary"); 
return(0););
for (i=2,n-1,
for (j=2,n-1,
for (k=j+1,n,
if(areprop(v[1]-v[i],v[j]-v[k])==1, print(Delta "is bad: " v[1]","v[i]","v[j]","v[k]" have rational cross-ratio");return(0););
);););
print(Delta " is good!");
);
return(1);
}




/*--------------------------------------
can difference of roots of of f divided by that of g be rational?
input: Q-irreducible monic polynomials f,g in Z[x] of degree at least 2; we assume that a root of f generates a Galois extension of Q, and the roots of g belong to this extension;
output: 0 if some difference of two roots of f divided by some difference of two roots of g   is (non-trivially) rational; 1 if none.  
--------------------------------------*/



ratrootdifquot(f,g)=
{
my(u,v,w,m,n,K,h);




m=poldegree(f);
n=poldegree(g);




if (m<2, print("deg f<2"); return(0););
if (n<2, print("deg g<2"); return(0););

if(polisirreducible(f)==0, print("f is reducible"); return(0););
if(polisirreducible(g)==0, print("g is reducible"); return(0););

h=subst(f,x,y); 

K=nfinit(h);
u=nfgaloisconj(K);
if(length(u)<m,
print(f "not Galois"); 
return(0););

if (f==g, return(iscrossratiorational(f));); 

/* from now on f\ne g*/ 

v=nfroots(K,g);

if(length(v)<n,
print("internal error"); 
return(0););

/*this must never happen but who knows */
  

w=vector(n,i,lift(v[i]));

/*extracting the first composant of the polmod type */ 

 


for (i=2,m,
for (j=1,n-1,
for (k=j+1,n,
if(areprop(u[1]-u[i],w[j]-w[k])==1, print(u[1]","u[i]","w[j]","w[k]" have rational cross-ratio");return(0););
);););
return(1);
}


/*----------------------------------------------------------- 
verification for Section 8.3 
input: list of pairs of discriminants 
output: 1 if OK, 0 if not
---------------------------------------------------------------*/ 


sec83(v)=
{
my(n,u,Deltax,Deltay,f,g);

n=length(v);

for(i=1,n,

u=v[i];

Deltax=u[1];
Deltay=u[2];

print(Deltax","Deltay);

f=polclass(Deltax);
g=polclass(Deltay);


if(ratrootdifquot(f,g)==0, 
print(Deltax","Deltay" bad");
return(0););

);

return(1);

}



/*---------------------------------------------------
list for Section 8.3
---------------------------------------------------*/

list93= [[-96,-192],   [-96, -288],   [-120, -160] ,   [-120, - 280], [-120, -760],   [-160, -280],   [-160, - 760],    [-180, -240],   [-192,-288],   [-195, -520],   [-195, - 715],   [-280, -760],  [-340, - 595],   [-480, - 960],    [-520, - 715]]


