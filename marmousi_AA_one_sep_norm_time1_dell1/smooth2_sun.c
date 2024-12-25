/************************************************************/
/*                   smooth the velocity	                */
/*                extracted from SU by sxd                  */
/*               author: sxd   date:2006.1                  */
/************************************************************/
#include "math.h"
#include "string.h"
#include <stdlib.h>
#include <stdio.h>

#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#define PI (3.141592653589793)
#define INITIAL_TSU 999999
#define NULL	((void *)0)
#define EPS FLT_MIN

/***************************************************************************/
  //void *ealloc1 (size_t n1, size_t size);
//void **ealloc2 (size_t n1, size_t n2, size_t size);
//void ***ealloc3 (size_t n1, size_t n2, size_t n3, size_t size);
void *alloc1 (size_t n1, size_t size);
void **alloc2 (size_t n1, size_t n2, size_t size);
void ***alloc3 (size_t n1, size_t n2, size_t n3, size_t size);
int *alloc1int (size_t n1);
float *alloc1float (size_t n1);
float **alloc2float (size_t n1, size_t n2);
float ***alloc3float (size_t n1, size_t n2, size_t n3);

float *ealloc1float(size_t n1);
float **ealloc2float(size_t n1, size_t n2);
float ***ealloc3float(size_t n1, size_t n2, size_t n3);

void free1float (float *p);
void free2float (float **p);
void free3float (float ***p);

void tripd (float *d, float *e, float *b, int n);

/****************************************************************/
/****************************************************************/
main()
{
	int n1;		/* number of points in x1 (fast) dimension */
	int n2;		/* number of points in x1 (fast) dimension */
	int nmax;	/* max of n1 and n2 */
	int ix, iz;	/* counters */
	int win[4];	/* 1d array defining the corners of smoothing window */
	float **v0;	/* array of input velocities */
	float **v;	/* array of output velocities */
	float **w;	/* intermediate array */
	float *errz;	/* array of error estimates as a function of x1 */
	float *d, *e;	/* input arrays for subroutine tripd */
	float *f;	/* intermediate array */
	float r1;	/* smoothing parameter for x1 direction */
	float r2;	/* smoothing parameter for x2 direction */
	float rw;	/* smoothing parameter for window */
	float err0;	/* error variable */
	float vrms;	/* rms velocity */
	FILE *infp;	/* input file pointer */
	FILE *outfp;	/* output file pointer */
	FILE *errorfp;		/* error file pointer */
	char *infile,*outfile,*errorfile;	/* name of error file */

/************input parameters ******************/
	infile="aoxian100_100.dat";
	outfile="vpsmooth.dat";
	errorfile="error.dat";
	n1=100;  /* NZ */
	n2=100;   /* NX */
	r1=5;     /*  1<=r1<=20 */
	r2=5;     /*  1<=r2<=20 */
    rw =0.2;

	win[0] = 0;
	win[1] = n1;
	win[2] = 0;
	win[3] = n2;


/***********************************************/
    infp= fopen(infile,"rb");
    outfp= fopen(outfile,"wb");
	errorfp= fopen(errorfile,"wb");
	/* scale the smoothing parameter */
	r1 = r1*r1*0.25;
	r2 = r2*r2*0.25;

	/* allocate space */
	nmax = (n1 < n2)? n2:n1;

	v = alloc2float(n1,n2);
	v0 = alloc2float(n1,n2);
	w = alloc2float(n1,n2);
	errz = alloc1float(nmax);
	d = alloc1float(nmax);
	e = alloc1float(nmax);
	f = alloc1float(nmax);

	/* read velocities */
	fread(v[0],sizeof(float),n2*n1,infp);

	/* save the original velocity */
        for(ix=0; ix<n2; ++ix)
	 	for(iz=0; iz<n1; ++iz)
			v0[ix][iz]=v[ix][iz];

	/* get parameters for window function */

	rw = rw*rw*0.25;

	/* define the window function */
	for(ix=0; ix<n2; ++ix)
	 	for(iz=0; iz<n1; ++iz)
			w[ix][iz] = 0;
	for(ix=win[2]; ix<win[3]; ++ix)
	 	for(iz=win[0]; iz<win[1]; ++iz)
			w[ix][iz] = 1;

	if(win[0]>0 || win[1]<n1 || win[2]>0 || win[3]<n2){
	/*	smooth the window function */
         	for(iz=0; iz<n1; ++iz){
	 		for(ix=0; ix<n2; ++ix){
				d[ix] = 1.0+2.0*rw;
				e[ix] = -rw;
				f[ix] = w[ix][iz];
			}
        		d[0] -= rw;
         		d[n2-1] -= rw;
         		tripd(d,e,f,n2);
	 		for(ix=0; ix<n2; ++ix)
				w[ix][iz] = f[ix];
		}
         	for(ix=0; ix<n2; ++ix){
	 		for(iz=0; iz<n1; ++iz){
				d[iz] = 1.0+2.0*rw;
				e[iz] = -rw;
				f[iz] = w[ix][iz];
		}
        		d[0] -= rw;
         		d[n1-1] -= rw;
         		tripd(d,e,f,n1);
	 		for(iz=0; iz<n1; ++iz)
				w[ix][iz] = f[iz];
		}
	}

	/*      solving for the smoothing velocity */
        for(iz=0; iz<n1; ++iz){
	 	for(ix=0; ix<n2-1; ++ix){
			d[ix] = 1.0+r2*(w[ix][iz]+w[ix+1][iz]);
			e[ix] = -r2*w[ix+1][iz];
			f[ix] = v[ix][iz];
		}
        	d[0] -= r2*w[0][iz];
         	d[n2-1] = 1.0+r2*w[n2-1][iz];
		f[n2-1] = v[n2-1][iz];
         	tripd(d,e,f,n2);
	 	for(ix=0; ix<n2; ++ix)
			v[ix][iz] = f[ix];
	}
         for(ix=0; ix<n2; ++ix){
	 	for(iz=0; iz<n1-2; ++iz){
			d[iz] = 1.0+r1*(w[ix][iz+1]+w[ix][iz+2]);
			e[iz] = -r1*w[ix][iz+2];
			f[iz] = v[ix][iz+1];
		}
		f[0] += r1*w[ix][1]*v[ix][0];
         	d[n1-2] = 1.0+r1*w[ix][n1-1];
		f[n1-2] = v[ix][n1-1];
         	tripd(d,e,f,n1-1);
	 	for(iz=0; iz<n1-1; ++iz)
			v[ix][iz+1] = f[iz];
	}
	/* write smoothed data */
	fwrite(v[0],sizeof(float),n1*n2,outfp);

	/* if the user specifies the name of a an error file*/




		/*calculate the RMS error of velocity */
/*
	err0 = 0.;
		vrms = 0;
      	 	 for(iz=0; iz<n1; ++iz){
		     for(ix=0; ix<n2; ++ix){
			  err0 += (v0[ix][iz]-v[ix][iz])*(v0[ix][iz]-v[ix][iz]);
			  vrms += v0[ix][iz]*v0[ix][iz];
		     }
			if (vrms != 0.0)
		 		errz[iz] = sqrt(err0/vrms);
			else
		 		errz[iz] = sqrt(err0/EPS);
		}
		fwrite(errz,sizeof(float),n1,errorfp);
		fclose(errorfp);
*/
   exit(0);

}
/***********************************************************************/
void *alloc1 (size_t n1, size_t size)
{
	void *p;

	if ((p=malloc(n1*size))==NULL)
		return NULL;
	return p;
}
void **alloc2 (size_t n1, size_t n2, size_t size)
{
	size_t i2;
	void **p;

	if ((p=(void**)malloc(n2*sizeof(void*)))==NULL)
		return NULL;
	if ((p[0]=(void*)malloc(n2*n1*size))==NULL) {
		free(p);
		return NULL;
	}
	for (i2=0; i2<n2; i2++)
		p[i2] = (char*)p[0]+size*n1*i2;
	return p;
}
void ***alloc3 (size_t n1, size_t n2, size_t n3, size_t size)
{
	size_t i3,i2;
	void ***p;

	if ((p=(void***)malloc(n3*sizeof(void**)))==NULL)
		return NULL;
	if ((p[0]=(void**)malloc(n3*n2*sizeof(void*)))==NULL) {
		free(p);
		return NULL;
	}
	if ((p[0][0]=(void*)malloc(n3*n2*n1*size))==NULL) {
		free(p[0]);
		free(p);
		return NULL;
	}

	for (i3=0; i3<n3; i3++) {
		p[i3] = p[0]+n2*i3;
		for (i2=0; i2<n2; i2++)
			p[i3][i2] = (char*)p[0][0]+size*n1*(i2+n2*i3);
	}
	return p;
}
int *alloc1int(size_t n1)
{
	return (int*)alloc1(n1,sizeof(int));
}
float *alloc1float(size_t n1)
{
	return (float*)alloc1(n1,sizeof(float));
}
float **alloc2float(size_t n1, size_t n2)
{
	return (float**)alloc2(n1,n2,sizeof(float));
}
float ***alloc3float(size_t n1, size_t n2, size_t n3)
{
	return (float***)alloc3(n1,n2,n3,sizeof(float));
}

float *ealloc1float(size_t n1)
{
	float *p;

	if (NULL == (p=alloc1float(n1)))
	{printf("%s: malloc failed");getchar();}
	return p;
}
float **ealloc2float(size_t n1, size_t n2)
{
	float **p;

	if (NULL == (p=alloc2float(n1, n2)))
		{printf("%s: malloc failed");getchar();}
	return p;
}
float ***ealloc3float(size_t n1, size_t n2, size_t n3)
{
	float ***p;

	if (NULL == (p=alloc3float(n1, n2, n3)))
	{printf("%s: malloc failed");getchar();}
	return p;
}

void free1 (void *p)
{
	free(p);
}
void free2 (void **p)
{
	free(p[0]);
	free(p);
}
void free3 (void ***p)
{
	free(p[0][0]);
	free(p[0]);
	free(p);
}
void free1float(float *p)
{
	free1(p);
}
void free2float(float **p)
{
	free2((void**)p);
}
void free3float(float ***p)
{
	free3((void***)p);
}

void tripd(float *d, float *e, float *b, int n)
/*****************************************************************************
Given an n-by-n symmetric, tridiagonal, positive definite matrix A and
 n-vector b, following algorithm overwrites b with the solution to Ax = b

*****************************************************************************
Input:
d[]	the diagonal of A
e[]	the superdiagonal of A
b[]	the rhs of Ax=b

Output:
b[]	b[] is overwritten with the solution to Ax=b

*****************************************************************************
Notes:

Given an n-by-n symmetric, tridiagonal, positive definite matrix A and
 n-vector b, following algorithm overwrites b with the solution to Ax = b

*****************************************************************************
Author: Zhenyue Liu, Colorado School of Mines, 1993.
*****************************************************************************/
{
	int k;
	float temp;

	/* decomposition */
	for(k=1; k<n; ++k){
           temp = e[k-1];
           e[k-1] = temp/d[k-1];
           d[k] -= temp*e[k-1];
	}

	/* substitution	*/
        for(k=1; k<n; ++k)  b[k] -= e[k-1]*b[k-1];

        b[n-1] /=d[n-1];
        for(k=n-1; k>0; --k)  b[k-1] = b[k-1]/d[k-1] - e[k-1]*b[k];

 }
