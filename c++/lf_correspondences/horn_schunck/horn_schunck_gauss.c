/*****************************************************************************/
/*                                                                           */
/*                   Copyright 08/2006 by Dr. Andres Bruhn,                  */
/*     Faculty of Mathematics and Computer Science, Saarland University,     */
/*                           Saarbruecken, Germany.                          */
/*                                                                           */
/*****************************************************************************/


#ifndef OF_HORN_SCHUNCK_INCLUDED
#define OF_HORN_SCHUNCK_INCLUDED

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc_mem_linear.c"
#include "alloc_mem_linear_mult.c"
#include "funct_lib.c"

#define TYPE_GRADIENT 0
#define TYPE_GRAY 1
#define TYPE_HESSIAN 2

/* ------------------------------------------------------------------------- */

void backward_registration
(
                           /**************************************************/
    float **f1,            /* in  : 1st image                                */
    float **f2,            /* in  : 2nd image                                */
    float **f2_bw,         /* out : 2nd image (motion compensated)           */
    float **u,             /* in     : x-component of displacement field     */
    float **v,             /* in     : y-component of displacement field     */
    int   nx,              /* in     : size in x-direction                   */
    int   ny,              /* in     : size in y-direction                   */
    int   bx,              /* in     : boundary size in x-direction          */
    int   by,              /* in     : boundary size in y-direction          */
    float hx,              /* in     : grid spacing in x-direction           */
    float hy               /* in     : grid spacing in y-direction           */
                           /**************************************************/
    )

/* creates warped version of image f2 by means of bilinear interpolation */

{
                           /**************************************************/
int   i,j;                 /* loop variables                                 */
int   ii,jj;               /* pixel coordinates                              */
float ii_fp,jj_fp;         /* subpixel coordinates                           */
float delta_i,delta_j;     /* subpixel displacement                          */
float hx_1,hy_1;           /* time saver                                     */
                           /**************************************************/

/* compute time savers */
hx_1=1.0/hx;
hy_1=1.0/hy;

/* set boundaries zero */
set_bounds_2d(f2,nx,ny,bx,by,(float)0.0);

 for (i=bx; i<nx+bx; i++)
     for (j=by; j<ny+by; j++)
     {
	 /* compute subpixel location */
	 ii_fp=i+(u[i][j]*hx_1);
	 jj_fp=j+(v[i][j]*hy_1);

	 /* if the required image information is out of bounds */
	 if ((ii_fp<bx)||(jj_fp<by)||(ii_fp>(nx+bx-1))||(jj_fp>(ny+by-1)))
	   {
	       /* set flow field zero*/
	       f2_bw[i][j]=0;

	   }
	 /* if required image information is available */
	 else
	   {
	       /* compute integer pixel coordinates of upper left pixel */

	       /* ----- TODO: fill in your code here ---- */

		   /* --------------------------------------- */
		ii = (int) floor(ii_fp);
		jj = (int) floor(jj_fp);


	       /* compute subpixel displacement */

	       /* ----- TODO: fill in your code here ---- */

		   /* --------------------------------------- */
	          delta_i = ii_fp-(float)ii;
		        delta_j = jj_fp-(float)jj;

	       /* perform bilinear interpolation */

	       /* ----- TODO: fill in your code here ---  */
		       f2_bw[i][j] = (1-delta_i)*(1-delta_j)*f2[ii][jj] + (delta_i)*(1-delta_j)*f2[ii+1][jj] + (1-delta_i)*delta_j*f2[ii][jj+1]  + delta_i*delta_j*f2[ii+1][jj+1];

		   /* --------------------------------------- */
	   }
     }
}


/********* */
float flow_driven_smooth_term(float **u, float **v,int i, int j, float hx, float hy, float epsilon){
  float ux, uy, vx,vy;
  float s=0;
  ux=(u[i+1][j]-u[i-1][j])/(2*hx);
  uy=(u[i][j+1]-u[i][j-1])/(2*hy);
  vx=(v[i+1][j]-v[i-1][j])/(2*hx);
  vy=(v[i][j+1]-v[i][j-1])/(2*hy);
  s= ux*ux+uy*uy+vx*vx+vy*vy;
  return 1.0/(2*sqrt(s+epsilon));

}


/* ------------------------------------------------------------------------- */

void horn_schunck_gaussseidel_direction_flow_driven
(
                        /*****************************************************/
float **J_11,           /* in     : entry 11 of the motion tensor            */
float **J_22,           /* in     : entry 22 of the motion tensor            */
float **J_33,           /* in     : entry 33 of the motion tensor            */
float **J_12,           /* in     : entry 12 of the motion tensor            */
float **J_13,           /* in     : entry 13 of the motion tensor            */
float **J_23,           /* in     : entry 23 of the motion tensor            */
float **u,              /* in+out : x-component of displacement field        */
float **v,              /* in+out : y-component of displacement field        */
int   nx,               /* in     : size in x-direction                      */
int   ny,               /* in     : size in y-direction                      */
int   bx,               /* in     : boundary size in x-direction             */
int   by,               /* in     : boundary size in y-direction             */
float hx,               /* in     : grid spacing in x-direction              */
float hy,               /* in     : grid spacing in y-direction              */
float s,                /* in     : norminated view direction x coordinates             */
float t,                /* in     : norminated view direction y coordinates             */
float beta,             /* in     : smoothness weight for direction */
float alpha             /* in     : smoothness weight                        */
                        /*****************************************************/
)

/*
 Computes one Gauss-Seidel iteration
*/

{
                        /*****************************************************/
int   i,j;              /* loop variables                                    */
float hx_2,hy_2;        /* time saver variables                              */
float xp,xm,yp,ym;      /* neighbourhood weights                             */
float sum;              /* central weight                                    */
float J_11_s, J_22_s, J_12_s; /* new motion tensor with direction */
float sxp,sxm,syp,sym; /* smoothness for neighbourhood */
float sxy; /* smoothness at current position */
                        /*****************************************************/


/* define time saver variables */
hx_2=alpha/(hx*hx);
hy_2=alpha/(hy*hy);

float epsilon = 0.0001;

/* set boundaries zero */
set_bounds_2d(u,nx,ny,bx,by,0.0);
set_bounds_2d(v,nx,ny,bx,by,0.0);


for(i=bx+1;i<nx+bx-1;i++)
    for(j=by+1;j<ny+by-1;j++)
    {
	/* compute weights */


  sxp = flow_driven_smooth_term(u,v,i+1,j,hx,hy,epsilon);
  sxm = flow_driven_smooth_term(u,v,i-1,j,hx,hy,epsilon);
  syp = flow_driven_smooth_term(u,v,i,j+1,hx,hy,epsilon);
  sym = flow_driven_smooth_term(u,v,i,j-1,hx,hy,epsilon);
  sxy = flow_driven_smooth_term(u,v,i,j,hx,hy,epsilon);

  xp =  (i<nx+bx-1) * hx_2 *(sxp+sxy)/2.0;
	xm =  (i>bx)      * hx_2 *(sxm+sxy)/2.0;
	yp =  (j<ny+by-1) * hy_2 *(syp+sxy)/2.0;
	ym =  (j>by)      * hy_2 *(sym+sxy)/2.0;

  J_11_s = J_11[i][j]+beta*t*t;
  J_22_s = J_22[i][j]+beta*s*s;
  J_12_s = J_12[i][j]-beta*s*t;

	sum = - (xp + xm + yp + ym);

	/* perform iteration */
	u[i][j]= (   J_13[i][j]
                   + J_12_s*v[i  ][j  ]
                         - xm * u[i-1][j  ]
                         - ym * u[i  ][j-1]
                         - yp * u[i  ][j+1]
		         - xp * u[i+1][j  ])
	        /(sum-J_11_s);

	v[i][j]= (   J_23[i][j]
                   + J_12_s*u[i  ][j  ]
                         - xm * v[i-1][j  ]
                         - ym * v[i  ][j-1]
                         - yp * v[i  ][j+1]
		         - xp * v[i+1][j  ])
	        /(sum-J_22_s);

    }

}




/* ------------------------------------------------------------------------- */

void horn_schunck_gaussseidel_direction
(
                        /*****************************************************/
float **J_11,           /* in     : entry 11 of the motion tensor            */
float **J_22,           /* in     : entry 22 of the motion tensor            */
float **J_33,           /* in     : entry 33 of the motion tensor            */
float **J_12,           /* in     : entry 12 of the motion tensor            */
float **J_13,           /* in     : entry 13 of the motion tensor            */
float **J_23,           /* in     : entry 23 of the motion tensor            */
float **u,              /* in+out : x-component of displacement field        */
float **v,              /* in+out : y-component of displacement field        */
int   nx,               /* in     : size in x-direction                      */
int   ny,               /* in     : size in y-direction                      */
int   bx,               /* in     : boundary size in x-direction             */
int   by,               /* in     : boundary size in y-direction             */
float hx,               /* in     : grid spacing in x-direction              */
float hy,               /* in     : grid spacing in y-direction              */
float s,                /* in     : norminated view direction x coordinates             */
float t,                /* in     : norminated view direction y coordinates             */
float beta,             /* in     : smoothness weight for direction */
float alpha             /* in     : smoothness weight                        */
                        /*****************************************************/
)

/*
 Computes one Gauss-Seidel iteration
*/

{
                        /*****************************************************/
int   i,j;              /* loop variables                                    */
float hx_2,hy_2;        /* time saver variables                              */
float xp,xm,yp,ym;      /* neighbourhood weights                             */
float sum;              /* central weight                                    */
float J_11_s, J_22_s, J_12_s; /* new motion tensor with direction */
                        /*****************************************************/


/* define time saver variables */
hx_2=alpha/(hx*hx);
hy_2=alpha/(hy*hy);


/* set boundaries zero */
set_bounds_2d(u,nx,ny,bx,by,0.0);
set_bounds_2d(v,nx,ny,bx,by,0.0);


for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {
	/* compute weights */
	xp =  (i<nx+bx-1) * hx_2;
	xm =  (i>bx)      * hx_2;
	yp =  (j<ny+by-1) * hy_2;
	ym =  (j>by)      * hy_2;
  J_11_s = J_11[i][j]+beta*t*t;
  J_22_s = J_22[i][j]+beta*s*s;
  J_12_s = J_12[i][j]-beta*s*t;

	sum = - (xp + xm + yp + ym);

	/* perform iteration */
	u[i][j]= (   J_13[i][j]
                   + J_12_s*v[i  ][j  ]
                         - xm * u[i-1][j  ]
                         - ym * u[i  ][j-1]
                         - yp * u[i  ][j+1]
		         - xp * u[i+1][j  ])
	        /(sum-J_11_s);

	v[i][j]= (   J_23[i][j]
                   + J_12_s*u[i  ][j  ]
                         - xm * v[i-1][j  ]
                         - ym * v[i  ][j-1]
                         - yp * v[i  ][j+1]
		         - xp * v[i+1][j  ])
	        /(sum-J_22_s);

    }

}


/* ------------------------------------------------------------------------- */

void horn_schunck_gaussseidel
(
                        /*****************************************************/
float **J_11,           /* in     : entry 11 of the motion tensor            */
float **J_22,           /* in     : entry 22 of the motion tensor            */
float **J_33,           /* in     : entry 33 of the motion tensor            */
float **J_12,           /* in     : entry 12 of the motion tensor            */
float **J_13,           /* in     : entry 13 of the motion tensor            */
float **J_23,           /* in     : entry 23 of the motion tensor            */
float **u,              /* in+out : x-component of displacement field        */
float **v,              /* in+out : y-component of displacement field        */
int   nx,               /* in     : size in x-direction                      */
int   ny,               /* in     : size in y-direction                      */
int   bx,               /* in     : boundary size in x-direction             */
int   by,               /* in     : boundary size in y-direction             */
float hx,               /* in     : grid spacing in x-direction              */
float hy,               /* in     : grid spacing in y-direction              */
float alpha             /* in     : smoothness weight                        */
                        /*****************************************************/
)

/*
 Computes one Gauss-Seidel iteration
*/

{
                        /*****************************************************/
int   i,j;              /* loop variables                                    */
float hx_2,hy_2;        /* time saver variables                              */
float xp,xm,yp,ym;      /* neighbourhood weights                             */
float sum;              /* central weight                                    */
                        /*****************************************************/


/* define time saver variables */
hx_2=alpha/(hx*hx);
hy_2=alpha/(hy*hy);


/* set boundaries zero */
set_bounds_2d(u,nx,ny,bx,by,0.0);
set_bounds_2d(v,nx,ny,bx,by,0.0);


for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {
	/* compute weights */
	xp =  (i<nx+bx-1) * hx_2;
	xm =  (i>bx)      * hx_2;
	yp =  (j<ny+by-1) * hy_2;
	ym =  (j>by)      * hy_2;

	sum = - (xp + xm + yp + ym);

	/* perform iteration */
	u[i][j]= (   J_13[i][j]
                   + J_12[i][j]*v[i  ][j  ]
                         - xm * u[i-1][j  ]
                         - ym * u[i  ][j-1]
                         - yp * u[i  ][j+1]
		         - xp * u[i+1][j  ])
	        /(sum-J_11[i][j]);

	v[i][j]= (   J_23[i][j]
                   + J_12[i][j]*u[i  ][j  ]
                         - xm * v[i-1][j  ]
                         - ym * v[i  ][j-1]
                         - yp * v[i  ][j+1]
		         - xp * v[i+1][j  ])
	        /(sum-J_22[i][j]);

    }

}


/* ------------------------------------------------------------------------- */

void compute_motion_tensor
(
                        /*****************************************************/
float **f1,             /* in     : 1st image                                */
float **f2,             /* in     : 2nd image                                */
int   nx,               /* in     : size in x-direction                      */
int   ny,               /* in     : size in y-direction                      */
int   bx,               /* in     : boundary size in x-direction             */
int   by,               /* in     : boundary size in y-direction             */
float hx,               /* in     : grid spacing in x-direction              */
float hy,               /* in     : grid spacing in y-direction              */
float **J_11,           /* out    : entry 11 of the motion tensor            */
float **J_22,           /* out    : entry 22 of the motion tensor            */
float **J_33,           /* out    : entry 33 of the motion tensor            */
float **J_12,           /* out    : entry 12 of the motion tensor            */
float **J_13,           /* out    : entry 13 of the motion tensor            */
float **J_23            /* out    : entry 23 of the motion tensor            */
                        /*****************************************************/
)

/*
 Computes the motion tensor entries from the given image pair
*/

{
                        /*****************************************************/
int     i,j;            /* loop variables                                    */
float   **fx;           /* first order image derivatives                     */
float   **fy;           /* first order image derivatives                     */
float   **ft;           /* first order image derivatives                     */
float   **fxx;          /* second order image derivatives                    */
float   **fxy;          /* second order image derivatives                    */
float   **fyy;          /* second order image derivatives                    */
float   **fxt;          /* second order image derivatives                    */
float   **fyt;          /* second order image derivatives                    */
float   hx_1,hy_1;      /* time saver variables                              */
                        /*****************************************************/


/* allocate memory */
ALLOC_MATRIX(8, nx+2*bx, ny+2*by,
	     &fx,
	     &fy,
	     &ft,
	     &fxx,
	     &fxy,
	     &fyy,
	     &fxt,
	     &fyt);

/* define time saver variables */
hx_1=1.0/(2.0*hx);
hy_1=1.0/(2.0*hy);

/* mirror boundaries */
mirror_bounds_2d(f1,nx,ny,bx,by);
mirror_bounds_2d(f2,nx,ny,bx,by);

/* compute first oder derivatives */
for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {
	fx[i][j] = 0.5*(f1[i+1][j]-f1[i-1][j]+f2[i+1][j]-f2[i-1][j])*hx_1;
	fy[i][j] = 0.5*(f1[i][j+1]-f1[i][j-1]+f2[i][j+1]-f2[i][j-1])*hy_1;
	ft[i][j] = (f2[i][j]-f1[i][j]);
    }

/* mirror boundaries */
mirror_bounds_2d(fx,nx,ny,bx,by);
mirror_bounds_2d(fy,nx,ny,bx,by);
mirror_bounds_2d(ft,nx,ny,bx,by);

/* compute second order derivatives */
for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {

	/* ----- TODO: fill in your code here ---- */
	fxx[i][j] = (fx[i+1][j]-fx[i-1][j])*hx_1;
	fyy[i][j] = (fy[i][j+1]-fy[i][j-1])*hy_1;
	fxy[i][j] = (fx[i][j+1]-fx[i][j-1])*hy_1;
	fxt[i][j] = hx_1*(ft[i+1][j] - ft[i-1][j]);//(-f1[i+1][j]+f1[i-1][j]+f2[i+1][j]-f2[i-1][j])*hx_1;
	fyt[i][j] = hy_1*(ft[i][j+1]-ft[i][j-1]);//(-f1[i][j+1]+f1[i][j-1]+f2[i][j+1]-f2[i][j-1])*hy_1;;

	/* --------------------------------------- */

    }


/* compute motion tensor entries entries */
for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {
	/* set up motion tensor for the gradient constancy assumption */


	/* ----- TODO: fill in your code here ---- */
	J_11[i][j] = fxx[i][j]*fxx[i][j] + fxy[i][j]*fxy[i][j];
	J_22[i][j] = fyy[i][j]*fyy[i][j] + fxy[i][j]*fxy[i][j];
	J_33[i][j] = fxt[i][j]*fxt[i][j] + fyt[i][j]*fyt[i][j];
	J_12[i][j] = fxx[i][j]*fxy[i][j] + fxy[i][j]*fyy[i][j];
	J_13[i][j] = fxt[i][j]*fxx[i][j] + fyt[i][j]*fxy[i][j];
	J_23[i][j] = fxt[i][j]*fxy[i][j] + fyt[i][j]*fyy[i][j];

	/* --------------------------------------- */

	/* --- TODO: then, remove tensor for the  grey value constancy
	             assumption                                          --- */

	/* grey value constancy assumption */

//	J_11[i][j] = fx[i][j] * fx[i][j];
//	J_22[i][j] = fy[i][j] * fy[i][j];
//	J_33[i][j] = ft[i][j] * ft[i][j];
//	J_12[i][j] = fx[i][j] * fy[i][j];
//	J_13[i][j] = fx[i][j] * ft[i][j];
//	J_23[i][j] = fy[i][j] * ft[i][j];

    /* --------------------------------------- */

    }

/* free memory */
FREE_MATRIX(8, nx+2*bx, ny+2*by,
	    fx,
	    fy,
	    ft,
	    fxx,
	    fxy,
	    fyy,
	    fxt,
	    fyt);
}


/* ------------------------------------------------------------------------- */


void HORN_SCHUNCK_VIEW
(
                        /*****************************************************/
float **f1,             /* in     : 1st image                                */
float **f2,             /* in     : 2nd image                                */
float **u,              /* out    : x-component of displacement field        */
float **v,              /* out    : y-component of displacement field        */
int   nx,               /* in     : size in x-direction                      */
int   ny,               /* in     : size in y-direction                      */
int   bx,               /* in     : boundary size in x-direction             */
int   by,               /* in     : boundary size in y-direction             */
float hx,               /* in     : grid spacing in x-direction              */
float hy,               /* in     : grid spacing in y-direction              */
float s,                /* in     : norminated view direction x coordinates             */
float t,                /* in     : norminated view direction y coordinates             */
float m_beta,             /* in     : smoothness weight for direction */
float m_alpha,          /* in     : smoothness weight                        */
int   n_iter            /* in     : number of iterations                     */
                        /*****************************************************/
)

/* computes optic flow with Horn/Schunck */

{

                        /*****************************************************/
int   i,j;              /* loop variables                                    */
float **J_11;           /* entry 11 of the motion tensor                     */
float **J_22;           /* entry 22 of the motion tensor                     */
float **J_33;           /* entry 33 of the motion tensor                     */
float **J_12;           /* entry 12 of the motion tensor                     */
float **J_13;           /* entry 13 of the motion tensor                     */
float **J_23;           /* entry 23 of the motion tensor                     */
                        /*****************************************************/



/* ---- alloc memory ---- */
ALLOC_MATRIX (6, nx+2*bx,  ny+2*by,
	      &J_11,
	      &J_22,
	      &J_33,
	      &J_12,
	      &J_13,
	      &J_23);


/* ---- initialise displacement field with zero ---- */
for (i=bx; i<bx+nx; i++)
    for (j=by; j<by+ny; j++)
    {
	u[i][j]=0;
	v[i][j]=0;
    }


/* ---- compute motion tensor ---- */
compute_motion_tensor(f1, f2, nx, ny, bx, by, hx, hy,
		      J_11, J_22, J_33, J_12, J_13, J_23);



/* ---- perform Gauss-Seidel iterations ---- */
for(i=1;i<=n_iter;i++)
 {
     horn_schunck_gaussseidel_direction_flow_driven(J_11, J_22, J_33, J_12, J_13, J_23,
			 u, v, nx, ny, bx, by, hx, hy,s,t,m_beta, m_alpha);
 }



/* ---- free memory ---- */
FREE_MATRIX (6, nx+2*bx,  ny+2*by,
	     J_11,
	     J_22,
	     J_33,
	     J_12,
	     J_13,
	     J_23);


}

void HORN_SCHUNCK
(
                        /*****************************************************/
float **f1,             /* in     : 1st image                                */
float **f2,             /* in     : 2nd image                                */
float **u,              /* out    : x-component of displacement field        */
float **v,              /* out    : y-component of displacement field        */
int   nx,               /* in     : size in x-direction                      */
int   ny,               /* in     : size in y-direction                      */
int   bx,               /* in     : boundary size in x-direction             */
int   by,               /* in     : boundary size in y-direction             */
float hx,               /* in     : grid spacing in x-direction              */
float hy,               /* in     : grid spacing in y-direction              */
float m_alpha,          /* in     : smoothness weight                        */
int   n_iter            /* in     : number of iterations                     */
                        /*****************************************************/
)

/* computes optic flow with Horn/Schunck */

{

                        /*****************************************************/
int   i,j;              /* loop variables                                    */
float **J_11;           /* entry 11 of the motion tensor                     */
float **J_22;           /* entry 22 of the motion tensor                     */
float **J_33;           /* entry 33 of the motion tensor                     */
float **J_12;           /* entry 12 of the motion tensor                     */
float **J_13;           /* entry 13 of the motion tensor                     */
float **J_23;           /* entry 23 of the motion tensor                     */
                        /*****************************************************/



/* ---- alloc memory ---- */
ALLOC_MATRIX (6, nx+2*bx,  ny+2*by,
	      &J_11,
	      &J_22,
	      &J_33,
	      &J_12,
	      &J_13,
	      &J_23);


/* ---- initialise displacement field with zero ---- */
for (i=bx; i<bx+nx; i++)
    for (j=by; j<by+ny; j++)
    {
	u[i][j]=0;
	v[i][j]=0;
    }


/* ---- compute motion tensor ---- */
compute_motion_tensor(f1, f2, nx, ny, bx, by, hx, hy,
		      J_11, J_22, J_33, J_12, J_13, J_23);



/* ---- perform Gauss-Seidel iterations ---- */
for(i=1;i<=n_iter;i++)
 {
     horn_schunck_gaussseidel(J_11, J_22, J_33, J_12, J_13, J_23,
			 u, v, nx, ny, bx, by, hx, hy, m_alpha);
 }



/* ---- free memory ---- */
FREE_MATRIX (6, nx+2*bx,  ny+2*by,
	     J_11,
	     J_22,
	     J_33,
	     J_12,
	     J_13,
	     J_23);


}
/* ------------------------------------------------------------------------- */

#endif
