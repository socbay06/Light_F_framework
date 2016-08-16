/*****************************************************************************/
/*                                                                           */
/*                   Copyright 08/2006 by Dr. Andres Bruhn,                  */
/*     Faculty of Mathematics and Computer Science, Saarland University,     */
/*                           Saarbruecken, Germany.                          */
/*                                                                           */
/*****************************************************************************/

#ifndef MG_TRANS_LIB_INCLUDED
#define MG_TRANS_LIB_INCLUDED

#include <math.h>
#include "alloc_mem_linear.c"
#include "conv_lib.c"


/* ------------------------------------------------------------------------- */


void resample_1d

(
                /*************************************************************/
    float *u,   /* in     : input vector, size 1..n                          */
    int   n,    /* in     : size of input vector                             */
    int   m,    /* in     : size of output vector                            */
    float *v    /* out    : output vector, size 1..m                         */
                /*************************************************************/
    )

/* Area-based resampling: Transforms a 1D image u of size n into an image v  */
/* of size m by integration over piecewise constant functions. Conservative. */

{
                         /****************************************************/
int     i, k;            /* loop variables                                   */
float   hu, hv;          /* grid sizes                                       */
float   uleft, uright;   /* boundaries                                       */
float   vleft, vright;   /* boundaries                                       */
float   fac;             /* normalization factor                             */
                         /****************************************************/

/*****************************************************************************/
/* (1/2) Special cases of area-based resampling are computed efficiently     */
/*****************************************************************************/

/* fast interpolation for output images of even size */ 
if (m==2*n)
 {  
     /* one cell is devided in two cells with equal value */
     for (i=1; i<=n; i++)
     {       
     v[i*2-1]=u[i];
     v[i*2  ]=u[i];
     }
     return;
 }


/* fast restriction for input images of even size */
if (2*m==n)
 {  
     /* two celss are melted to a larger cell with averaged value */
     for (i=1; i<=m; i++)
     {
     v[i]=0.5*(u[i*2-1]+u[i*2]);
     }
     
     return;
     }


/*****************************************************************************/
/* (2/2) Remaining cases require more complex algorithm                      */
/*****************************************************************************/


/* initializations */
                            /*************************************************/
hu    = 1.0 / (float)n;     /* grid size of u                                */
hv    = 1.0 / (float)m;     /* grid size of v                                */
uleft = 0.0;                /* left interval boundary of u                   */
vleft = 0.0;                /* left interval boundary of v                   */
k     = 1;                  /* index for u                                   */
fac   = hu / hv;            /* for normalization                             */
                            /*************************************************/


/*---- loop ----*/

for (i=1; i<=m; i++)
    /* calculate v[i] by integrating the piecewise constant function u */
    {
    /* calculate right interval boundaries */
    uright = uleft + hu;
    vright = vleft + hv;

    if (uright > vright)
       /* since uleft <= vleft, the entire v-cell i is in the u-cell k */
       v[i] = u[k];
    else
       {
       /* consider fraction alpha of the u-cell k in v-cell i */
       v[i] = (uright - vleft) * n * u[k++];
     
       /* update */
       uright = uright + hu;

       /* consider entire u-cells inside v-cell i */
       while (uright <= vright)
             /* u-cell k lies entirely in v-cell i; sum up */
             {
             v[i] = v[i] + u[k++];
             uright = uright + hu;
             }

       /* consider fraction beta of the u-cell k in v-cell i */
       v[i] = v[i] + (1.0 - (uright - vright) * n) * u[k];

       /* normalization */
       v[i] = v[i] * fac;
       } /* else */

    /* update */
    uleft = uright - hu;
    vleft = vright;
    /* now it holds: uleft <= vleft */
    }  /* for i */


return;

}


/*--------------------------------------------------------------------------*/



void resample_2d_x

(
                /*************************************************************/
float  **u,     /* in   : input image                                        */
int    nx,      /* in   : x dimension of input image                         */
int    ny,      /* in   : y dimension of input image                         */
int    bx,      /* in   : size of border in x-direction                      */
int    by,      /* in   : size of border in y-direction                      */
float  **u_out, /* out  : output image                                       */
int    mx       /* in   : x dimension of output image                        */
                /*************************************************************/
)

/* resample a 2-D image in x-direction using area-based resampling */

{
                         /****************************************************/
int    i, j;             /* loop variables                                   */
float  *uhelp, *vhelp;   /* auxiliary vectors                                */
                         /****************************************************/


/* allocate memory */
ALLOC_VECTOR (1, nx+2, &uhelp);
ALLOC_VECTOR (1, mx+2, &vhelp);

/* resample image linewise in x-direction */
for (j=by; j<ny+by; j++)
 {       
     /* initialise left boundary of 1-D array with zero */
     uhelp[0]=u[bx][j];

     /* copy current line in this 1-D array */
     for (i=bx; i<nx+bx; i++) 
	 uhelp[i-bx+1] = u[i][j];
     
    /* initialise right boundary of this 1-D array with zero */
    uhelp[nx+1]=u[nx+bx-1][j];

    /* resample this 1-D array */
    resample_1d (uhelp, nx, mx, vhelp);

    /* copy resmapled array in corresponding output line */
    for (i=bx; i<mx+bx; i++)
	u_out[i][j] = vhelp[i-bx+1];
    }

/* free memory */
FREE_VECTOR (1, nx+2, uhelp);
FREE_VECTOR (1, mx+2, vhelp);

return;
} 

/*---------------------------------------------------------------------------*/


void resample_2d_y

(
                /*************************************************************/
float  **u,     /* in   : input image                                        */
int    nx,      /* in   : x dimension of input image                         */
int    ny,      /* in   : y dimension of input image                         */
int    bx,      /* in   : size of border in x-direction                      */
int    by,      /* in   : size of border in y-direction                      */
float  **u_out, /* out  : output image                                       */
int    my       /* in   : y dimension of output image                        */
                /*************************************************************/
)

/* resample a 2-D image in y-direction using area-based resampling */

{
                         /****************************************************/
int    i, j;             /* loop variables                                   */
float  *uhelp, *vhelp;   /* auxiliary vectors                                */
                         /****************************************************/


/* allocate memory */
ALLOC_VECTOR (1, ny+2, &uhelp);
ALLOC_VECTOR (1, my+2, &vhelp);

/* resample image columnwise in y-direction */
for (i=bx; i<nx+bx; i++)
 {       
     /* initialsie left boundary of 1-D array with zero */
     uhelp[0]=u[i][by];

     /* copy current column in this 1-D array */
     for (j=by; j<ny+by; j++) 
	 uhelp[j-by+1] = u[i][j];

     /* initialise right boundary of this 1-D array with zero */
     uhelp[ny+1]=u[i][ny+by-1];
     
     /* resample this 1-D array */
     resample_1d (uhelp, ny, my, vhelp);

     /* copy resmapled array in corresponding output column */
     for (j=by; j<my+by; j++) 
	 u_out[i][j] = vhelp[j-by+1];
     
 }

/* free memory */
FREE_VECTOR (1, ny+2, uhelp);
FREE_VECTOR (1, my+2, vhelp);

return;
}


/*--------------------------------------------------------------------------*/



void resample_2d

(
                /*************************************************************/
float  **u,     /* in   : input image                                        */
int    nx,      /* in   : x dimension of input image                         */
int    ny,      /* in   : y dimension of input image                         */
int    bx,      /* in   : size of border in x-direction                      */
int    by,      /* in   : size of border in y-direction                      */
float  **u_out, /* out  : output image                                       */
int    mx,      /* in   : x dimension of output image                        */
int    my,      /* in   : y dimension of output image                        */
float  **tmp    /* tmp  : temporary array of size (nx * my)                  */
                /*************************************************************/
)

/* resample a 2-D image using area-based resampling */

{

    /* if interpolation */
    if (my>=ny)
    {
	/* resample first in x-direction */
	resample_2d_x(u, nx, ny, bx, by, tmp, mx);
	/* resample then in y-direction */
	resample_2d_y(tmp, mx, ny, bx, by, u_out, my);
    }
    /* if restriction */
    else 
    {
	/*  resample first in y-direction */
	resample_2d_y(u, nx, ny, bx, by, tmp, my);
        /*  resample tehn in x-direction */
	resample_2d_x(tmp, nx, my, bx, by, u_out, mx);
    }

return;

}


/*--------------------------------------------------------------------------*/

#endif
