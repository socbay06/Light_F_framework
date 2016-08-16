/*****************************************************************************/
/*                                                                           */
/*                   Copyright 08/2006 by Dr. Andres Bruhn,                  */
/*     Faculty of Mathematics and Computer Science, Saarland University,     */
/*                           Saarbruecken, Germany.                          */
/*                                                                           */
/*****************************************************************************/


#ifndef CONV_LIB_INCLUDED
#define CONV_LIB_INCLUDED

#include <math.h>
#include <stdlib.h>
#include "alloc_mem_linear.c"
#include "matrix_lib.c"

/* ------------------------------------------------------------------------- */


void gauss_mask_odd

(
                       /******************************************************/
    float   sigma,     /* in     : standard deviation of Gaussian            */
    int     nm,        /* in     : data dimension in mask direktion          */
    float   hm,        /* in     : pixel size in mask direction              */
    float   precision, /* in     : cutoff at precision * sigma               */
    int     *masksize, /* out    : size of gaussian mask                     */
    float   **mask     /* out    : the mask itself                           */
                       /*          memory is allocated IN the routine !      */
                       /******************************************************/
)

/* computes odd gaussian convolution mask */

{
                       /******************************************************/
    int    i;          /* loop variables                                     */
    float  sum;        /* for summing up                                     */
                       /******************************************************/


/* calculate size of convolution mask n                                     */
/* please note that due to the symmetry of the Gaussian convolution mask    */
/* only n+1 instead of 2*n+1 value have to be stored                        */
*masksize = (int)(precision * sigma / hm) + 1;
if ((*masksize) > nm)
   {
       printf("\n gauss_mask_odd : mask too large");
       printf("\n nm   : %d",nm);
       printf("\n hm   : %f",hm);
       printf("\n size : %d",*masksize);
       printf("\n prec : %f\n\n",precision);
       exit(0);
   }

/* allocate storage for convolution vector (n+1 values) */
ALLOC_VECTOR (1, (*masksize)+1, mask);


/* calculate entries of convolution vector */
 for (i=0; i<=(*masksize); i++)
   {
     (*mask)[i] = 1 / (sigma * sqrt(2.0 * 3.1415926)) 
	 * exp (- (i * i * hm * hm) / (2.0 * sigma * sigma));
   }


/* normalise convolution vector to sum 1 */
 sum = (*mask)[0];
 for (i=1; i<=(*masksize); i++)
     sum = sum + 2.0 * (*mask)[i];
 for (i=0; i<=(*masksize); i++)
     (*mask)[i] = (*mask)[i] / sum;

 
 return;
}



/* ------------------------------------------------------------------------- */


void conv_2d_x_sym_odd_opt

(
                       /******************************************************/
    int     masksize,  /* in     : size of convolution mask                  */
    float   *mask,     /* in     : convolution mask                          */
    int     nx,        /* in     : data dimension in x direction             */
    int     ny,        /* in     : data dimension in y direction             */
    int     bx,        /* in     : boundary in x direction                   */
    int     by,        /* in     : boundary in y direction                   */
    float   **f,       /* in     : original data                             */
    float   **v        /* out    : processed data                            */
                       /******************************************************/
)

/* convolution in x-direction with odd symmetric convolution mask */
/* since the values are stored in y-direction in the cache, a loop unrolling */
/* scheme is applied */

{

                            /*************************************************/
int    i, j, k,p;           /* loop variables                                */
float  *sum;                /* for summing up                                */
float  **help;              /* array of rows with suitable boundary size     */
int    tmp1,tmp2,tmp3,tmp4; /* time saver                                    */
int    UNROLL;              /* number of rows that are convolved in parallel */
int    inner_loop_max;      /* number of loops for parrallel computation     */
int    inner_loop_rest;     /* number of remaining rows                      */
                            /*************************************************/


/* set number of rows convolved in parallel */
UNROLL=32; 

/* allocate storage for that many rows */
ALLOC_MATRIX (1, nx+masksize+masksize,UNROLL, &help);

/* allocate storagy for that many results */
ALLOC_VECTOR (1, UNROLL, &sum);

/* compute number of loops required if the desired number of rows is */
/* convolved in parallel */
inner_loop_max=ny/UNROLL;

/* compute number of remaining rows that have to be convolved thereafter */
inner_loop_rest=ny-(inner_loop_max*UNROLL);

/* time saver indices */ 
tmp1=masksize-1;
tmp2=tmp1+nx;
tmp3=tmp2+1;
tmp4=tmp1-(bx-1);

/*****************************************************************************/
/* (1/2) As long as the desired number of rows can be convolved in parallel  */
/*       use loop unrolling scheme that respects cache direction             */
/*****************************************************************************/

 for (j=0; j<inner_loop_max; j++)
 {
     /* copy rows in vector array */
     for (i=bx; i<nx+bx; i++)
	 for (k=0; k<UNROLL; k++)
	 {
	     help[i+tmp4][k] = f[i][j*UNROLL+k+by];
	 }
     
     /* mirror boundaries of each of these rows */
     for (p=1; p<=masksize; p++)
	 for (k=0; k<UNROLL; k++)
	 {
	     help[masksize-p][k] = help[tmp1+p][k];
	     help[tmp2    +p][k] = help[tmp3-p][k];
	 }
     
     /* convolution step for each of these rows */
     for (i=masksize; i<=tmp2; i++)
     {
	 /* convolve different rows in parallel */
	 for (k=0; k<UNROLL; k++)
	 {	    
	     sum[k] = mask[0] * help[i][k];
	 }
	 
	 for (p=1; p<=masksize; p++)
	     for (k=0; k<UNROLL; k++)
	     {
		 sum[k] += mask[p] * (help[i+p][k] + help[i-p][k]);
	     } 
	 
	 /* write back results in parallel */
	 for (k=0; k<UNROLL; k++)	
	 {	
	     v[i-tmp4][j*UNROLL+k+by] = sum[k];     
	 }
     }
 } /* for j */

/*****************************************************************************/
/* (2/2) Convolve the remaining number of rows in parallel using the same    */
/*       loop unrolling scheme                                               */
/*****************************************************************************/

 if (inner_loop_rest>0)
 {
     /* copy rows in vector array */
     for (i=bx; i<nx+bx; i++)
	 for (k=0; k<inner_loop_rest; k++)
	 {
	     help[i+tmp4][k] = f[i][j*UNROLL+k+by];
	 }
     
     /* mirror boundaries for each of these rows */
     for (p=1; p<=masksize; p++)
	 for (k=0; k<inner_loop_rest; k++)
	 {
	     help[masksize-p][k]      = help[tmp1+p][k];
	     help[tmp2+p][k] = help[tmp3-p][k];
	 }

     /* convolution step for each of these rows */
     for (i=masksize; i<=tmp2; i++)
     {
	 /* convolve different rows in parallel */
	 for (k=0; k<inner_loop_rest; k++)
	 {	    
	     sum[k] = mask[0] * help[i][k];
	 }

	 for (p=1; p<=masksize; p++)
	     for (k=0; k<inner_loop_rest; k++)
	     {
		 sum[k] += mask[p] * (help[i+p][k] + help[i-p][k]);
	     } 

	  /* write back results in parallel */
	 for (k=0; k<inner_loop_rest; k++)
	 {	    
	     v[i-tmp4][j*UNROLL+k+by] = sum[k];
	 }
     }
 }

/* disallocate storage for the rows */
FREE_MATRIX (1, nx+masksize+masksize,UNROLL, help);

/* disallocate storage for the results */
FREE_VECTOR (1, UNROLL, sum);
 
return;
}




/* ------------------------------------------------------------------------- */


void conv_2d_y_sym_odd_opt


(
                       /******************************************************/
    int     masksize,  /* in     : size of convolution mask                  */
    float   *mask,     /* in     : convolution mask                          */
    int     nx,        /* in     : data dimension in x direction             */
    int     ny,        /* in     : data dimension in y direction             */
    int     bx,        /* in     : boundary in x direction                   */
    int     by,        /* in     : boundary in y direction                   */
    float   **f,       /* in     : original data                             */
    float   **v        /* out    : processed data                            */
                       /******************************************************/
)

/* convolution in y-direction with odd symmetric convolution mask */

{
                            /*************************************************/
int    i, j, k,p;           /* loop variables                                */
float  sum;                 /* for summing up                                */
float  *help;               /* array for one column with suitable boundary   */
int    tmp1,tmp2,tmp3,tmp4; /* time saver                                    */
                            /*************************************************/


/* allocate storage for a single row */
ALLOC_VECTOR (1, ny+masksize+masksize, &help);

/* time saver indices */
tmp1=masksize-1;
tmp2=tmp1+ny;
tmp3=tmp2+1;
tmp4=tmp1-(by-1);


 /* for each column */
 for (i=bx; i<nx+bx; i++)
 {
     /* copy current column in column vector */
     for (j=by; j<ny+by; j++)
	 help[j+tmp4] = f[i][j];

     /* mirror boundaries of the column vector */
     for (p=1; p<=masksize; p++)
     {
	 help[masksize-p]      = help[tmp1+p];
	 help[tmp2+p] = help[tmp3-p];
     }

     /* convolution step */
     for (j=masksize; j<=tmp2; j++)
     {
	 /* calculate convolution */
	 sum = mask[0] * help[j];
	 for (p=1; p<=masksize; p++)
	     sum += mask[p] * (help[j+p] + help[j-p]);
	 /* write back result of current column */
	 v[i][j-tmp4] = sum;
     }
 } /* for i */

/* disallocate storage for a single row */
FREE_VECTOR (1, ny+masksize+masksize, help);

return;

}


/* ------------------------------------------------------------------------- */


void conv_2d_x_asym_odd_opt

(
                       /******************************************************/
    int     masksize,  /* in     : size of convolution mask                  */
    float   *mask,     /* in     : convolution mask                          */
    int     nx,        /* in     : data dimension in x direction             */
    int     ny,        /* in     : data dimension in y direction             */
    int     bx,        /* in     : boundary in x direction                   */
    int     by,        /* in     : boundary in y direction                   */
    float   **f,       /* in     : original data                             */
    float   **v        /* out    : processed data                            */
                       /******************************************************/
)

/* convolution in x-direction with odd symmetric convolution mask */
/* since the values are stored in y-direction in the cache, a loop unrolling */
/* scheme is applied */

{

                            /*************************************************/
int    i, j, k,p;           /* loop variables                                */
float  *sum;                /* for summing up                                */
float  **help;              /* array of rows with suitable boundary size     */
int    tmp1,tmp2,tmp3,tmp4; /* time saver                                    */
int    UNROLL;              /* number of rows that are convolved in parallel */
int    inner_loop_max;      /* number of loops for parrallel computation     */
int    inner_loop_rest;     /* number of remaining rows                      */
                            /*************************************************/


/* set number of rows convolved in parallel */
UNROLL=32; 

/* allocate storage for that many rows */
ALLOC_MATRIX (1, nx+masksize+masksize,UNROLL, &help);

/* allocate storagy for that many results */
ALLOC_VECTOR (1, UNROLL, &sum);

/* compute number of loops required if the desired number of rows is */
/* convolved in parallel */
inner_loop_max=ny/UNROLL;

/* compute number of remaining rows that have to be convolved thereafter */
inner_loop_rest=ny-(inner_loop_max*UNROLL);

/* time saver indices */ 
tmp1=masksize-1;
tmp2=tmp1+nx;
tmp3=tmp2+1;
tmp4=tmp1-(bx-1);


/*****************************************************************************/
/* (1/2) As long as the desired number of rows can be convolved in parallel  */
/*       use loop unrolling scheme that respects cache direction             */
/*****************************************************************************/

 for (j=0; j<inner_loop_max; j++)
 {
     /* copy row in vector array */
     for (i=bx; i<nx+bx; i++)
	 for (k=0; k<UNROLL; k++)
	 {
	     help[i+tmp4][k] = f[i][j*UNROLL+k+by];
	 }
     
     /* mirror boundaries for each of these rows */
     for (p=1; p<=masksize; p++)
	 for (k=0; k<UNROLL; k++)
	 {
	     help[masksize-p][k]      = help[tmp1+p][k];
	     help[tmp2+p][k] = help[tmp3-p][k];
	 }

     /* convolution step for each of these rows */
     for (i=masksize; i<=tmp2; i++)
     {
	 /* convolve different rows in parallel */
	 for (k=0; k<UNROLL; k++)
	 {	   
	     sum[k] = mask[0] * help[i][k];
	 }

	 for (p=1; p<=masksize; p++)
	     for (k=0; k<UNROLL; k++)
	     {
		 sum[k] += mask[p] * (help[i+p][k] - help[i-p][k]);
	     } 

	 /* write back results in parallel */
	 for (k=0; k<UNROLL; k++)
	 {	    
	     v[i-tmp4][j*UNROLL+k+by] = sum[k];
	 }
     }
 } /* for j */

/*****************************************************************************/
/* (2/2) Convolve the remaining number of rows in parallel using the same    */
/*       loop unrolling scheme                                               */
/*****************************************************************************/

 if (inner_loop_rest>0)
 {
     /* copy rows in vector array */
     for (i=bx; i<nx+bx; i++)
	 for (k=0; k<inner_loop_rest; k++)
	 {
	     help[i+tmp4][k] = f[i][j*UNROLL+k+by];
	 }
     
     /* mirror boundaries for each of these rows */
     for (p=1; p<=masksize; p++)
	 for (k=0; k<inner_loop_rest; k++)
	 {
	     help[masksize-p][k]      = help[tmp1+p][k];
	     help[tmp2+p][k] = help[tmp3-p][k];
	 }

     /* convolution step for each of these rows */
     for (i=masksize; i<=tmp2; i++)
     {
	 /* convolve different rows in parallel */
	 for (k=0; k<inner_loop_rest; k++)
	 {	     
	     sum[k] = mask[0] * help[i][k];
	 }

	 for (p=1; p<=masksize; p++)
	     for (k=0; k<inner_loop_rest; k++)
	     {
		 sum[k] += mask[p] * (help[i+p][k] - help[i-p][k]);
	     } 

	 /* write back results in parallel */
	 for (k=0; k<inner_loop_rest; k++)
	 {	     
	     v[i-tmp4][j*UNROLL+k+by] = sum[k];
	 }
     }
 }

/* disallocate storage for the rows */
FREE_MATRIX (1, nx+masksize+masksize,UNROLL, help);

/* disallocate storage for the results */
FREE_VECTOR (1, UNROLL, sum);
 
return;
}


/* ------------------------------------------------------------------------- */


void conv_2d_y_asym_odd_opt

(
                       /******************************************************/
    int     masksize,  /* in     : size of convolution mask                  */
    float   *mask,     /* in     : convolution mask                          */
    int     nx,        /* in     : data dimension in x direction             */
    int     ny,        /* in     : data dimension in y direction             */
    int     bx,        /* in     : boundary in x direction                   */
    int     by,        /* in     : boundary in y direction                   */
    float   **f,       /* in     : original data                             */
    float   **v        /* out    : processed data                            */
                       /******************************************************/
)

/* convolution in y-direction with odd antisymmetric convolution mask */

{
                            /*************************************************/
int    i, j, k,p;           /* loop variables                                */
float  sum;                 /* for summing up                                */
float  *help;               /* array for one column with suitable boundary   */
int    tmp1,tmp2,tmp3,tmp4; /* time saver                                    */
                            /*************************************************/


/* allocate storage for a sigle row */
ALLOC_VECTOR (1, ny+masksize+masksize, &help);

/* time saver indices */
tmp1=masksize-1;
tmp2=tmp1+ny;
tmp3=tmp2+1;
tmp4=tmp1-(by-1);

 /* for each column */
 for (i=bx; i<nx+bx; i++)
 {
     /* copy current column in this column vector */
     for (j=by; j<ny+by; j++)
	 help[j+tmp4] = f[i][j];

     /* mirror boundary of the colum vector */
     for (p=1; p<=masksize; p++)
     {
	 help[masksize-p]      = help[tmp1+p];
	 help[tmp2+p] = help[tmp3-p];
     }
     
     /* convolution step for the column vector*/
     for (j=masksize; j<=tmp2; j++)
     {
	 /* calculate convolution */
	 sum = mask[0] * help[j];
	 for (p=1; p<=masksize; p++)
	     sum += mask[p] * (help[j+p] - help[j-p]);
	 /* write back result for the current colum */
	 v[i][j-tmp4] = sum;
     }
 } /* for i */  

/* disallocate storage for a single row */
FREE_VECTOR (1, ny+masksize+masksize, help);

return;

}

/* ------------------------------------------------------------------------- */

     
void presmooth_2d

(
                     /********************************************************/
float **f,           /* in  : original image                                 */
float **f_smoothed,  /* out : smoothed image                                 */
int  nx,             /* in  : size in x-direction                            */
int  ny,             /* in  : size in y-direction                            */
int  bx,             /* in  : boundary size in x-direction                   */
int  by,             /* in  : boundary size in y-direction                   */
float hx,            /* in  : grid size in x-direction                       */
float hy,            /* in  : grid size in y-direction                       */
float sigmax,        /* in  : std. dev. of Gaussian in y-direction           */
float sigmay         /* in  : std. dev. of Gaussian in z-direction           */
                     /********************************************************/
)

/* performs convolution of image f with Gaussian of standard deviation */
/* sigmax and sigmay in x- and y-direction, respectively */

{
   
                     /********************************************************/
float *maskx;        /* Gaussian convolution mask in x-direction             */
float *masky;        /* Gaussian convolution mask in y-direction             */
int   masksizex;     /* size of convolution mask in x-direction              */
int   masksizey;     /* size of convolution mask in y-direction              */
                     /********************************************************/

 /* compute Gaussian convolution mask in x- and y-direction */  
 if (sigmax>0) gauss_mask_odd(sigmax ,nx, hx ,(float)3.0,&masksizex,&maskx);
 if (sigmay>0) gauss_mask_odd(sigmay ,ny, hy ,(float)3.0,&masksizey,&masky);
 
 /* perform convolution with Gaussians in x- and y-direction */
 if (sigmax>0)  
     conv_2d_x_sym_odd_opt(masksizex,maskx,nx,ny,bx,by,f,f_smoothed);
 if (sigmax==0) 
     copy_matrix_2d(f,f_smoothed,nx,ny,bx,by);
 if (sigmay>0)  
     conv_2d_y_sym_odd_opt(masksizey,masky,nx,ny,bx,by,f_smoothed,f_smoothed);
 
 /* free memory allocated for Gaussian convolutions masks */
 if (sigmax>0) free(maskx);
 if (sigmay>0) free(masky);
}

/* ------------------------------------------------------------------------- */

#endif











