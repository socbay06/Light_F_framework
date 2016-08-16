#ifndef OF_LIB_INCLUDED
#define OF_LIB_INCLUDED

#include <math.h>
#include "alloc_mem_linear_mult.c"
#include "matrix_lib.c"


/* ------------------------------------------------------------------------- */

void calculate_errors_2d

(
                      /*******************************************************/
      float   **uref, /* in  : x component of reference flow field           */
      float   **vref, /* in  : y component of reference flow field           */
      float   **u,    /* in  : x component of computed flow field            */
      float   **v,    /* in  : y component of computed flow field            */
      int     nx,     /* in  : flow field  dimension in x direction          */
      int     ny,     /* in  : flow field  dimension in y direction          */
      int     bx,     /* in  : boundary size in x-direction                  */
      int     by,     /* in  : boundary size in y-direction                  */
      float   *aae,   /* out : average angular error                         */
      float   *al2e,  /* out : average l2 norm error                         */
      float   *ref_d, /* out : density of reference image                    */
      float   *cal_d  /* out : density of calculated image                   */
                      /*******************************************************/
)

//  calculates aeverage angular error and average l2 norm error

{
                    /*********************************************************/
float  sum1;        /* sum for average angular error                         */
float  sum2;        /* sum for average l2 norm error                         */
float  uc,vc;       /* reference ( (c)orrect ) flow field                    */
float  ue,ve;       /* computed  ( (e)stimated) flow field                   */
float  temp1,temp2; /* time saver                                            */
float  conv_const;  /* time saver                                            */
float  aux;         /* time saver                                            */
int    kcal;        /* percentage of pixels of the calculated flow field     */
                    /* where the optic flow exists                           */
int    kref;        /* percentage of pixels of the reference flow field      */
                    /* where the optic flow exists                           */
int    i,j;         /* loop variables                                        */
                    /*********************************************************/


/* initalise variables */
kcal = 0;
kref = 0;
sum1 = 0;
sum2 = 0;

/* set constant */
conv_const = (180.0 / 3.1415927) ;

 for (i=bx; i<nx+bx; i++)
     for (j=by; j<ny+by; j++)
     {
	 uc=uref[i][j];
	 vc=vref[i][j];
	 
	 /* if correct flow vector exists */
	 if ((uc != 100.0) || (vc != 100.0))
	 {
	     /* increase number of existing correct flow vectors */
	     kref = kref + 1;
	     
	     ue=u[i][j];
	     ve=v[i][j];

	     /* if estimated flow vector exist */
	     if ((ue != 100.0) || (ve != 100.0))
	     {  
                 /* increase number of existing estimated flow vectors */
		 kcal = kcal + 1;
		 
		 temp1=uc-ue;
		 temp2=vc-ve;
		 
		 /* compute l2 norm */
		 sum2=sum2 + sqrt(temp1*temp1 + temp2*temp2);
		 
		 /* compute angular error */
		 aux = (uc * ue + vc * ve + 1.0) 
		     / sqrt( (uc * uc + vc * vc + 1.0) 
			     * (ue * ue + ve * ve + 1.0) );
		 
		 if (aux > 1.0) aux = 1.0;
		 if (aux < - 1.0) aux = - 1.0;
		 
		 sum1=sum1 + conv_const * acos (aux);	    
	     }
	 }
     } 


/* calculate average norms */
*aae = sum1 / ((float) kcal);
*al2e = sum2 /((float) kcal);

/* calucalte densities */
*ref_d = (float)kref/(nx*ny);
*cal_d = (float)kcal/(nx*ny);

}


/* ------------------------------------------------------------------------- */

#endif
