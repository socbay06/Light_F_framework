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
#include "mg_trans_lib.c"


#define TYPE_GRADIENT 0
#define TYPE_GRAY 1
#define TYPE_HESSIAN 2

/* ------------------------------------------------------------------------- */

void backward_registration
(
        /**************************************************/
        float **f1,        /* in  : 1st image                                */
        float **f2,        /* in  : 2nd image                                */
        float **f2_bw,     /* out : 2nd image (motion compensated)           */
        float **u,         /* in     : x-component of displacement field     */
        float **v,         /* in     : y-component of displacement field     */
        int nx,            /* in     : size in x-direction                   */
        int ny,            /* in     : size in y-direction                   */
        int bx,            /* in     : boundary size in x-direction          */
        int by,            /* in     : boundary size in y-direction          */
        float hx,          /* in     : grid spacing in x-direction           */
        float hy           /* in     : grid spacing in y-direction           */
                           /**************************************************/
)

/* creates warped version of image f2 by means of bilinear interpolation */

{
        /**************************************************/
        int i,j;           /* loop variables                                 */
        int ii,jj;         /* pixel coordinates                              */
        float ii_fp,jj_fp; /* subpixel coordinates                           */
        float delta_i,delta_j; /* subpixel displacement                          */
        float hx_1,hy_1;   /* time saver                                     */
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
                                /* assume zero flow, i.e. set warped 2nd image to 1st image */
                                f2_bw[i][j]=f1[i][j];

                        }
                        /* if required image information is available */
                        else
                        {
                                /* compute integer index of upper left pixel */
                                ii=(int)floor(ii_fp);
                                jj=(int)floor(jj_fp);

                                /* compute subpixel displacement */
                                delta_i = ii_fp-(float)ii;
                                delta_j = jj_fp-(float)jj;

                                /* perform bilinear interpolation */
                                f2_bw[i][j]   = (1.0-delta_i)*(1.0-delta_j) * f2[ii  ][jj  ]
                                                +      delta_i *(1.0-delta_j) * f2[ii+1][jj  ]
                                                + (1.0-delta_i)*     delta_j  * f2[ii  ][jj+1]
                                                +      delta_i *     delta_j  * f2[ii+1][jj+1];
                        }
                }
}



/*--------------------------------------------------------------------------*/

void compute_max_warp_levels
(
        /*************************************************/
        int nx_orig,       /* in  : size in x-direction on finest grid      */
        int ny_orig,       /* in  : size in y-direction on finest grid      */
        float n_eta,       /* in  : warping reduction factor                */
        int   *max_levels  /* out : maximum number of levels                */
                           /*************************************************/
)

/* compute maximum number of warping levels for given image size and warping */
/* reduction factor */

{

        /*************************************************/
        int i;             /* level counter                                 */
        float nx,ny;       /* reduced dimensions                            */
        float nx_old,ny_old; /* reduced dimensions                            */
                             /*************************************************/

        nx_old=nx_orig;
        ny_old=ny_orig;

        /* Test different maximal number of warping levels i */
        for (i=1;; i++)
        {

                /* Compute corresponding resolution */
                nx=(int)ceil(nx_orig*pow(n_eta,i));
                ny=(int)ceil(ny_orig*pow(n_eta,i));

                /* if the resulting image size would be smaller than 4x4, we found the
                   maximal number of warp levels */
                if ((nx<4)||(ny<4)) break;

                nx_old=nx;
                ny_old=ny;
        }

        if((nx==1)||(ny==1)) i--;
        *max_levels=i;
}


/* ------------------------------------------------------------------------- */

void horn_schunck_warp_sor
(
        /*****************************************************/
        float **J_11,   /* in     : entry 11 of the motion tensor            */
        float **J_22,   /* in     : entry 22 of the motion tensor            */
        float **J_33,   /* in     : entry 33 of the motion tensor            */
        float **J_12,   /* in     : entry 12 of the motion tensor            */
        float **J_13,   /* in     : entry 13 of the motion tensor            */
        float **J_23,   /* in     : entry 23 of the motion tensor            */
        float **du,     /* in+out : x-component of flow increment            */
        float **dv,     /* in+out : y-component of flow increment            */
        float **u,      /* in     : x-component of flow field                */
        float **v,      /* in     : y-component of flow field                */
        int nx,         /* in     : size in x-direction                      */
        int ny,         /* in     : size in y-direction                      */
        int bx,         /* in     : boundary size in x-direction             */
        int by,         /* in     : boundary size in y-direction             */
        float hx,       /* in     : grid spacing in x-direction              */
        float hy,       /* in     : grid spacing in y-direction              */
        float alpha,    /* in     : smoothness weight                        */
        float omega     /* in     : SOR overrelaxation parameter             */
                        /*****************************************************/
)

/*
   Computes one SOR iteration
 */

{
        /*****************************************************/
        int i,j;        /* loop variables                                    */
        float hx_2,hy_2; /* time saver variables                              */
        float xp,xm,yp,ym; /* neighbourhood weights                             */
        float sum;      /* central weight                                    */
                        /*****************************************************/


/* define time saver variables */
        hx_2=alpha/(hx*hx);
        hy_2=alpha/(hy*hy);


/* set boundaries zero */
        set_bounds_2d(u,nx,ny,bx,by,0.0);
        set_bounds_2d(v,nx,ny,bx,by,0.0);
        set_bounds_2d(du,nx,ny,bx,by,0.0);
        set_bounds_2d(dv,nx,ny,bx,by,0.0);


        for(i=bx; i<nx+bx; i++)
                for(j=by; j<ny+by; j++)
                {
                        /* compute weights */
                        xp =  (i<nx+bx-1) * hx_2;
                        xm =  (i>bx)      * hx_2;
                        yp =  (j<ny+by-1) * hy_2;
                        ym =  (j>by)      * hy_2;

                        sum = xp + xm + yp + ym;

                        /* perform SOR iteration */

                        /* ----- TODO: fill in your code here ---- */

                        du[i][j] = (1-omega) * du[i][j] + omega *
                                   ( ( J_12[i][j] *dv[i][j] +J_13[i][j]
                                       - xm*( du[i-1][j] + u[i-1][j] )
                                       - xp*( du[i+1][j] + u[i+1][j] )
                                       - ym*( du[i][j-1] + u[i][j-1] )
                                       - yp*( du[i][j+1] + u[i][j+1] )
                                       + sum * ( u[i][j] ))
                                     /(-J_11[i][j] -sum ) );
                        dv[i][j] = (1-omega) * dv[i][j] + omega *
                                   ( ( J_12[i][j] *du[i][j] +J_23[i][j]
                                       - xm*( dv[i-1][j] + v[i-1][j] )
                                       - xp*( dv[i+1][j] + v[i+1][j] )
                                       - ym*( dv[i][j-1] + v[i][j-1] )
                                       - yp*( dv[i][j+1] + v[i][j+1] )
                                       + sum * ( v[i][j] ))
                                     /(-J_22[i][j] -sum ) );

                        /* --------------------------------------- */

                }
}

/* ------------------------------------------------------------------------- */

void compute_motion_tensor
(
        /*****************************************************/
        float **f1,     /* in     : 1st image                                */
        float **f2,     /* in     : 2nd image                                */
        int nx,         /* in     : size in x-direction                      */
        int ny,         /* in     : size in y-direction                      */
        int bx,         /* in     : boundary size in x-direction             */
        int by,         /* in     : boundary size in y-direction             */
        float hx,       /* in     : grid spacing in x-direction              */
        float hy,       /* in     : grid spacing in y-direction              */
        int type,       /* in     : which constancy assumption applied       */
        float **J_11,   /* out    : entry 11 of the motion tensor            */
        float **J_22,   /* out    : entry 22 of the motion tensor            */
        float **J_33,   /* out    : entry 33 of the motion tensor            */
        float **J_12,   /* out    : entry 12 of the motion tensor            */
        float **J_13,   /* out    : entry 13 of the motion tensor            */
        float **J_23    /* out    : entry 23 of the motion tensor            */
                        /*****************************************************/
)

/*
   Computes the motion tensor entries from a given image pair
 */

{
        /*****************************************************/
        int i,j;        /* loop variables                                    */
        float   **fx;   /* first order image derivatives                     */
        float   **fy;   /* first order image derivatives                     */
        float   **ft;   /* first order image derivatives                     */
        float   **fxx;  /* second order image derivatives                    */
        float   **fxy;  /* second order image derivatives                    */
        float   **fyy;  /* second order image derivatives                    */
        float   **fxt;  /* second order image derivatives                    */
        float   **fyt;  /* second order image derivatives                    */
        float hx_1,hy_1; /* time saver variables                              */
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
        for(i=bx; i<nx+bx; i++)
                for(j=by; j<ny+by; j++)
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
        for(i=bx; i<nx+bx; i++)
                for(j=by; j<ny+by; j++)
                {
                        fxx[i][j]=(fx[i+1][j  ]-fx[i-1][j  ])*hx_1;
                        fxy[i][j]=(fy[i+1][j  ]-fy[i-1][j  ])*hx_1;
                        fyy[i][j]=(fy[i  ][j+1]-fy[i  ][j-1])*hy_1;
                        fxt[i][j]=(ft[i+1][j  ]-ft[i-1][j  ])*hx_1;
                        fyt[i][j]=(ft[i  ][j+1]-ft[i  ][j-1])*hy_1;
                }


/* compute motion tensor entries entries */
        for(i=bx; i<nx+bx; i++)
                for(j=by; j<ny+by; j++)
                {
                        /* set up motion tensor */
                        if(type == TYPE_GRAY){

                        /* Brightness constancy */

                           J_11[i][j] = fx[i][j] * fx[i][j];
                           J_22[i][j] = fy[i][j] * fy[i][j];
                           J_33[i][j] = ft[i][j] * ft[i][j];
                           J_12[i][j] = fx[i][j] * fy[i][j];
                           J_13[i][j] = fx[i][j] * ft[i][j];
                           J_23[i][j] = fy[i][j] * ft[i][j];

                       }else{

                       /* Gradient constancy */
                       J_11[i][j] = fxx[i][j] * fxx[i][j] + fxy[i][j] * fxy[i][j];
                       J_22[i][j] = fxy[i][j] * fxy[i][j] + fyy[i][j] * fyy[i][j];
                       J_33[i][j] = fxt[i][j] * fxt[i][j] + fyt[i][j] * fyt[i][j];
                       J_12[i][j] = fxx[i][j] * fxy[i][j] + fxy[i][j] * fyy[i][j];
                       J_13[i][j] = fxx[i][j] * fxt[i][j] + fxy[i][j] * fyt[i][j];
                       J_23[i][j] = fxy[i][j] * fxt[i][j] + fyy[i][j] * fyt[i][j];
                     }
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

void HORN_SCHUNCK_WARP_LEVEL
(
        /*****************************************************/
        float **f1,     /* in     : 1st image                                */
        float **f2,     /* in     : 2nd image                                */
        float **du,     /* in+out : x-component of flow increment            */
        float **dv,     /* in+out : y-component of flow increment            */
        float **u,      /* in+out : x-component of flow field                */
        float **v,      /* in+out : y-component of flow field                */
        int nx,         /* in     : size in x-direction                      */
        int ny,         /* in     : size in y-direction                      */
        int bx,         /* in     : boundary size in x-direction             */
        int by,         /* in     : boundary size in y-direction             */
        float hx,       /* in     : grid spacing in x-direction              */
        float hy,       /* in     : grid spacing in y-direction              */
        int type,       /* in     : type of constancy assumption             */
        float m_alpha,  /* in     : smoothness weight                        */
        int n_iter,     /* in     : number of iterations                     */
        float n_omega   /* in     : SOR overrelaxation parameter             */
                        /*****************************************************/
)


{
        /*****************************************************/
        int i;          /* loop variable                                     */
        float  **J_11;  /* entry 11 of the motion tensor                     */
        float  **J_22;  /* entry 22 of the motion tensor                     */
        float  **J_33;  /* entry 33 of the motion tensor                     */
        float  **J_12;  /* entry 12 of the motion tensor                     */
        float  **J_13;  /* entry 13 of the motion tensor                     */
        float  **J_23;  /* entry 23 of the motion tensor                     */
                        /*****************************************************/



/* ---- alloc memory ---- */
        ALLOC_MATRIX (6, nx+2*bx,  ny+2*by,
                      &J_11,
                      &J_22,
                      &J_33,
                      &J_12,
                      &J_13,
                      &J_23);

/* ---- initialise displacement field with zero ----  */
        set_matrix_2d(du,nx+2*bx,ny+2*by,0,0,(float)0.0);
        set_matrix_2d(dv,nx+2*bx,ny+2*by,0,0,(float)0.0);


/* ---- compute motion tensor ---- */
        compute_motion_tensor(f1, f2, nx, ny, bx, by, hx, hy,type,
                              J_11, J_22, J_33, J_12, J_13, J_23);



/* ---- perform SOR iterations ---- */
        for(i=1; i<=n_iter; i++)
        {
                horn_schunck_warp_sor(J_11, J_22, J_33, J_12, J_13, J_23,
                                      du, dv, u, v, nx, ny, bx, by, hx, hy,
                                      m_alpha, n_omega);
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


void HORN_SCHUNCK_WARP
(
        /*****************************************************/
        float **f1_orig, /* in     : 1st image (original resolution)          */
        float **f2_orig, /* in     : 2nd image (original resolution)          */
        int nx_orig,    /* in     : size in x-direction (original resolution)*/
        int ny_orig,    /* in     : size in y-direction (original resoluiton)*/
        int type,       /* in     : type of constancy assumption             */
        float **f1_res, /* in+out : 1st image, resampled                     */
        float **f2_res, /* in+out : 2nd image, resampled                     */
        float **f2_res_warp, /* in+out : 2nd image, resampled  and warped         */
        float **du,     /* in+out : x-component of flow increment            */
        float **dv,     /* in+out : y-component of flow increment            */
        float **u,      /* in+out : x-component of flow field                */
        float **v,      /* in+out : y-component of flow field                */
        float **tmp,    /* in+out : temporary aray for resampling            */
        int nx_fine,    /* in     : size in x-direction (current resolution) */
        int ny_fine,    /* in     : size in y-direction (current resolution) */
        int bx,         /* in     : boundary size in x-direction             */
        int by,         /* in     : boundary size in y-direction             */
        float hx_fine,  /* in     : spacing in x-direction (current resol.)  */
        float hy_fine,  /* in     : spacing in y-direction (current resol.)  */
        float m_alpha,  /* in     : smoothness weight                        */
        int n_iter,     /* in     : number of iterations                     */
        float n_omega,  /* in     : SOR overrelaxation parameter             */
        float n_warp_eta, /* in     : warping reduction factor between levels  */
        int max_rec_depth, /* in     : maximum recursion depth                  */
        int rec_depth   /* in     : current recursion depth                  */
                        /*****************************************************/
)

/* implements warping for the Horn and Schunck method */

{

        /************************************************/
        int nx_coarse,ny_coarse; /* dimensions on previous coarser grid          */
        float hx_coarse,hy_coarse; /* grid sizes on previous coarser grid          */
                                   /************************************************/


/* compute dimensions and grid sizes for previous coarser grid */
        nx_coarse=(int)ceil(nx_orig*pow(n_warp_eta,rec_depth+1));
        ny_coarse=(int)ceil(ny_orig*pow(n_warp_eta,rec_depth+1));
        hx_coarse=(float)nx_orig/(float)nx_coarse;
        hy_coarse=(float)ny_orig/(float)ny_coarse;


/* start at coarsest level by recursively calling the routine */
        if (rec_depth < max_rec_depth)
        {
                HORN_SCHUNCK_WARP(f1_orig, f2_orig, nx_orig, ny_orig, type,
                                  f1_res, f2_res, f2_res_warp,
                                  du, dv, u, v, tmp,
                                  nx_coarse, ny_coarse, bx, by, hx_coarse, hy_coarse,
                                  m_alpha,
                                  n_iter, n_omega, n_warp_eta,
                                  max_rec_depth,rec_depth+1);
        }

/* ---- resample images ---------------------------------------------------- */

/* restrict original image pair to resolution of current level */
        resample_2d(f1_orig,nx_orig,ny_orig,bx,by,f1_res,nx_fine,ny_fine,tmp);
        resample_2d(f2_orig,nx_orig,ny_orig,bx,by,f2_res,nx_fine,ny_fine,tmp);


/* ---- get overall flow field from previous resolution level -------------- */

/* if on coarsest resolution */
        if(rec_depth==max_rec_depth)
        {
                /* set flow field zero */
                set_matrix_2d(u,nx_fine+2*bx,ny_fine+2*by,0,0,(float)0.0);
                set_matrix_2d(v,nx_fine+2*bx,ny_fine+2*by,0,0,(float)0.0);
        }
        /* if not on coarsest resolution */
        else
        {
                /* interpolate solution from previous coarser level */
                resample_2d(u,nx_coarse,ny_coarse,bx,by,u,nx_fine,ny_fine,tmp);
                resample_2d(v,nx_coarse,ny_coarse,bx,by,v,nx_fine,ny_fine,tmp);
        }


/* ---- set up difference problem at current resolution -------------------- */

/* warp second image by overall flow field from previous coarser resolution */
        backward_registration(f1_res,f2_res,f2_res_warp,u,v,
                              nx_fine,ny_fine,bx,by,hx_fine,hy_fine);



/* ---- solve difference problem at current resolution --------------------- */

/* solve difference problem at current resolution to obtain increment */
        HORN_SCHUNCK_WARP_LEVEL(f1_res, f2_res_warp, du, dv, u, v,
                                nx_fine, ny_fine, bx, by, hx_fine, hy_fine,type,
                                m_alpha, n_iter, n_omega);



/* ---- compute overall flow field at current resolution ------------------- */

/* sum up flow increment */
        add_matrix_2d(u,du,u,nx_fine,ny_fine,bx,by);
        add_matrix_2d(v,dv,v,nx_fine,ny_fine,bx,by);
}



/* ------------------------------------------------------------------------- */


void HORN_SCHUNCK_MAIN
(
        /*****************************************************/
        float **f1,     /* in     : 1st image                                */
        float **f2,     /* in     : 2nd image                                */
        float **u,      /* out    : x-component of displacement field        */
        float **v,      /* out    : y-component of displacement field        */
        int nx,         /* in     : size in x-direction                      */
        int ny,         /* in     : size in y-direction                      */
        int bx,         /* in     : boundary size in x-direction             */
        int by,         /* in     : boundary size in y-direction             */
        float hx,       /* in     : grid spacing in x-direction              */
        float hy,       /* in     : grid spacing in y-direction              */
        int type,       /* in     : type of constancy assumption             */
        float m_alpha,  /* in     : smoothness weight                        */
        int n_iter,     /* in     : number of iterations                     */
        float n_omega,  /* in     : SOR overrelaxation parameter             */
        float n_warp_eta, /* in     : warping reduction factor between levels  */
        int n_warp_levels /* in     : desired number of warping levels         */
                          /*****************************************************/
)

/* computes optic flow with Horn/Schunck + Warping */

{

        /*****************************************************/
        float **du;     /* x-component of flow increment                     */
        float **dv;     /* y-component of flow increment                     */
        float **f1_res; /* 1st image, resampled                              */
        float **f2_res; /* 2nd image, resampled                              */
        float **f2_res_warp; /* 2nd image, resampled  and warped                  */
        float **tmp;    /* temporary array for resampling                    */
        int max_rec_depth; /* maximum recursion depth (warping level -1)        */
        int n_warp_max_levels; /* maximum possible number of warping levels         */
                              /*****************************************************/


/* compute maximal number of warping levels given the downsampling factor eta
   and the image dimensions */
        compute_max_warp_levels(nx,ny,n_warp_eta,&n_warp_max_levels);

/* limit number of desired warping levels by number of possible levels  */
/* we have to substract one, since: n warping levels -> n-1 recursions  */
        max_rec_depth=minimum(n_warp_levels,n_warp_max_levels)-1;


/* ---- alloc memory ---- */
        ALLOC_MATRIX (6, nx+2*bx,  ny+2*by,
                      &du,
                      &dv,
                      &f1_res,
                      &f2_res,
                      &f2_res_warp,
                      &tmp);


/* call Horn/Schunck warping routine with desired number of levels */
        HORN_SCHUNCK_WARP(f1, f2, nx, ny,type,
                          f1_res, f2_res, f2_res_warp,
                          du, dv, u, v, tmp,
                          nx, ny, bx, by, hx, hy,
                          m_alpha,
                          n_iter, n_omega, n_warp_eta,
                          max_rec_depth,0);


/* ---- free memory ---- */
        FREE_MATRIX (6, nx+2*bx,  ny+2*by,
                     du,
                     dv,
                     f1_res,
                     f2_res,
                     f2_res_warp,
                     tmp);


}
/* ------------------------------------------------------------------------- */

#endif
