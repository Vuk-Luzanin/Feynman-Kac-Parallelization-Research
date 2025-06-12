#include <stdio.h>
#include <math.h>
#include <mpi.h>
/* sample size for each (x,y,z) initial point */
#define N 1000
/* number of (x,y,z) points */
#define NP 1000
#define twopi 6.2831853071795864770
/* these parameters define the 3-D ellipsoid */
#define a 3.0
#define b 2.0
#define c 1.0
#define v(p,q,r) 2.0*(p*p*am4+q*q*bm4+r*r*cm4)+am2+bm2+cm2
main(int argc,char **argv)
{
/* 
   Exercise 2.2 of Petersen and Arbenz, "Intro
   to Parallel Computing, Oxford Univ. Press,  2003.
   Problem is to solve PDE (Lap=Laplace operator in 3-D):

       (1/2) Lap u - v(x,y,z) u = 0,

   inside elliptic domain
 
       D = {(x,y,z)| x^2/a^2+y^2/b^2+z^2/c^2 <= 1}
   
   where u(x,y,z) = solution, u(boundary(D)) = 1.
   The potential v(x,y,z) is

       v = 2*(x^2/a^4+y^2/b^4+z^2/c^4) + 1/a^2 + 1/b^2 + 1/c^2.

   The analytic solution is u = exp(x^2/a^2+y^2/b^2+z^2/c^2 - 1).
   Our method is via Feynman-Kac Formula (see P&A Section 2.5,
   exercise 2.2). The idea is to start from any (x,y,z) in D, and
   compute (x+Wx(t),y+Wy(t),z+Wz(t)) where 3-D Brownian motion
   (Wx,Wy,Wz) is updated each step by sqrt(h)*(z1,z2,z3),
   each z1,z2,z3 are independent approximately Gaussian 
   random variables with zero mean and variance 1. Each 
   (x1(t),x2(t),x3(t)) is advanced until (x1,x2,x3) exits 
   the domain D. Upon its first exit from D, the sample 
   path (x1,x2,x3) is stopped and a new sample path at 
   (x,y,z) is started until N such paths are completed. 
   The Feynman-Kac formula gives the solution here as

       u(x,y,z) = (1/N) sum(i=1..N) Y(tau_i),

   where

       Y(tau) = exp( -int(s=0..tau) v(x1(s),x2(s),x3(s)) ds),

   and tau = first exit time for path (x1,x2,x3). Integration 
   procedure is a second order weak accurate method:

       X(t+h)  = [x1(t) + sqrt(h)*z1]
                 [x2(t) + sqrt(h)*z2]
                 [x3(t) + sqrt(h)*z3]

   again, the z1,z2,z3 are approximately normal univariate 
   Gaussians. An Euler predictor approximates Y at the end 
   of the step

       Y_e     = (1 - h*v(X(t)) * Y(t), 

   while a trapezoidal rule completes the step:

       Y(t+h)  = Y(t) - (h/2)*[v(X(t+h))*Y_e + v(X(t))*Y(t)].

                             wpp, 15 Dec. 2003, ETHZ

*/

    double x,y,z,x1,x2,x3,z1,z2,z3;
    double gerr,err,Y,Ye,Yt,soln,dx,dy,dz,vs,vh;
    double am2,bm2,cm2,am4,bm4,cm4,abcm2;
    double chk,fnm1,fnpm1,rt3h,sum,t3rd,ut,us,h=0.001; 
    int i,ierr,ip,npp,size,rank,master=0;
    int iprime[] = {331, 337, 347, 349, 353, 359, 367, 373, 379,
                    383, 389, 397, 401, 409, 419, 421, 431, 433};
    static double seed;
    double ggl(double*);
    MPI_Status stat;
/* end of declarations */

/* intialize MPI and get rank */
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    am2  = 1.0/(a*a); am4  = 1.0/(a*a*a*a);
    bm2  = 1.0/(b*b); bm4  = 1.0/(b*b*b*b);
    cm2  = 1.0/(c*c); cm4  = 1.0/(c*c*c*c);
    rt3h = sqrt(3.0*h);              /* step for B-increments */

/* loop over initial x,y,z points */

    err   = 0.0;
    t3rd  = 1.0/3.0;
    fnm1  = 1.0/((double) N);
    npp   = NP/size;      /* number of initial points/cpu */
    fnpm1 = 1.0/((double) npp);
    if(rank < 18){
/* if NCPU < 18, choose seed from iprime list */
        seed = (double) iprime[rank];
    } else {
/* if NCPU > 17, kludge "primes" off end of list, these
   are not all primes */
       ip   = (rank-16)/2;
       if((2*(ip/2))!=ip){
          seed = (double) (iprime[17] + ip*6);
       } else {
          seed = (double) (iprime[17] + ip*6 + 4);
       }
    }
    if(rank==master) gerr = 0.0;
    for(ip=0;ip<npp;ip++){

/* choose (x,y,z) point inside ellipse (x/a)^2+(y/b)^2+(z/c)^2 <= 1 */

       chk = 10000.0;
       while(chk > 1.0){
          x   = a*(2.0*ggl(&seed) - 1.0);   /* picks x in (-a,a) */
          y   = b*(2.0*ggl(&seed) - 1.0);   /* picks y in (-b,b) */
          z   = c*(2.0*ggl(&seed) - 1.0);   /* picks z in (-c,c) */
          chk = x*x*am2+y*y*bm2+z*z*cm2;    /* check if x,y,z in ellipse */
       }

/* compute exact solution for (x,y,z) */

       soln = exp(x*x*am2+y*y*bm2+z*z*cm2 - 1.0);
       Yt = 0.0;

/* do N paths */

       for(i=0;i<N;i++){
          x1  = x;
          x2  = y;
          x3  = z;
          Y   = 1.0;    /* Y = exp(-int(s=0..t) v(X)ds) */
          chk = 0.0;
          while(chk<1.0){
             ut = ggl(&seed);
             if(ut < t3rd){
                us = ggl(&seed) - 0.5;
                if(us < 0.0){dx = -rt3h;} else {dx = rt3h;}
             } else {
                dx = 0.0;
             }
             ut = ggl(&seed);
             if(ut < t3rd){
                us = ggl(&seed) - 0.5;
                if(us < 0.0){dy = -rt3h;} else {dy = rt3h;}
             } else {
                dy = 0.0;
             }
             ut = ggl(&seed);
             if(ut < t3rd){
                us = ggl(&seed) - 0.5;
                if(us < 0.0){dz = -rt3h;} else {dz = rt3h;}
             } else {
                dz = 0.0;
             }
             vs  = v(x1,x2,x3);
             x1 += dx;
             x2 += dy;
             x3 += dz;
             vh  = v(x1,x2,x3);
             Ye  = (1.0 - h*vs)*Y;
             Y   = Y - 0.5*h*(vh*Ye + vs*Y); 
             chk = x1*x1*am2+x2*x2*bm2+x3*x3*cm2;
          }
          Yt += Y;
       }
       Yt  *= fnm1; 
       err += (soln-Yt)*(soln-Yt);    /* square error of pt. ip */
    }
    err *= fnpm1;        /* mean square error of ip pts. */
    if(rank!=master){
       ierr  = MPI_Send(&err,1,MPI_DOUBLE,master,rank,MPI_COMM_WORLD);
       if(ierr!=0) printf(" MPI_Send from rank=%d fails\n",rank);
    } else {
       gerr += err;             /* include rank=0 computation */
       for(ip=1;ip<size;ip++){  /* loop over other CPU contributions */
          ierr = MPI_Recv(&err,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
                          MPI_COMM_WORLD,&stat);
          if(ierr!=0) printf(" MPI_Recv for ip=%d fails\n",ip);
          gerr += err;
       }
       gerr *= 1.0/((double) size);
       gerr  = sqrt(gerr);
       printf(" for NP=%d interior points, N=%d sample, NCPUs=%d\n",NP,N,size);
       printf(" npp=%d interior points/cpu, total=%d\n",npp,size*npp);
       printf("      h=%6.4f, MC solution RMS err = %e\n",h,gerr);
    }
    MPI_Finalize();
} 
#include <math.h>
double ggl(double *ds)
{
/* generate u(0,1) distributed random numbers.
   Seed ds must be saved between calls. ggl is
   essentially the same as the IMSL routine RNUM.

   W. Petersen and M. Troyer, 24 Oct. 2002, ETHZ */

   double t,d2=0.2147483647e10;
   t   = *ds;
   t   = fmod(0.16807e5*t,d2);
   *ds = t;
   return((t-1.0e0)/(d2-1.0e0));
}