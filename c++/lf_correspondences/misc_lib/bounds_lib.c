/*****************************************************************************/
/*                                                                           */
/*                   Copyright 08/2006 by Dr. Andres Bruhn                   */
/*     Faculty of Mathematics and Computer Science, Saarland University,     */
/*                           Saarbruecken, Germany.                          */
/*                                                                           */
/*****************************************************************************/


#ifndef BOUNDS_LIB_INCLUDED
#define BOUNDS_LIB_INCLUDED
   
/*---------------------------------------------------------------------------*/


void set_bounds_2d

(
                     /********************************************************/
    float **A,       /* image matrix                                         */
    int   nx,        /* size in x direction                                  */
    int   ny,        /* size in y direction                                  */
    int   bx,        /* boundary in x direction                              */
    int   by,        /* boundary in y direction                              */
    float a          /* set boundaries to a                                  */
                     /********************************************************/
)

/* set boundaries of the 2-D array A to value a */

{
                     /********************************************************/
int  i,j ;           /* loop variables                                       */
                     /********************************************************/

/* upper and lower boundary */
 for (i=0; i<nx+2*bx; i++)
     for (j=1; j<=by; j++)
     {
	 A[i][by   -j]=a;   
	 A[i][ny+by-1+j]=a;
     }
 
/* left and right boundary */
 for (i=1; i<=bx; i++)
     for (j=by; j<ny+by; j++)
     {
	 A[bx     -i][j]=a;   
	 A[nx+bx-1+i][j]=a;
     }
 
return;
}

/*---------------------------------------------------------------------------*/


void mirror_bounds_2d
(
                     /********************************************************/
    float **A,       /* image matrix                                         */
    int   nx,        /* size in x direction                                  */
    int   ny,        /* size in y direction                                  */
    int   bx,        /* boundary in x direction                              */
    int   by         /* boundary in y direction                              */
                     /********************************************************/
)

/* mirror boundaries of the 2-D array A */

{
                     /********************************************************/
int  i,j ;           /* loop variables                                       */
                     /********************************************************/

  
/* upper and lower boundary */
 for (i=bx; i<nx+bx; i++)
     for (j=1; j<=by; j++)
     {
	 A[i][by     -j]=A[i][by-1 +j];   
	 A[i][ny+by-1+j]=A[i][ny+by-j];
     }

/* left and right boundary */
 for (i=1; i<=bx; i++)
     for (j=0; j<ny+2*by; j++)
     {
	 A[bx     -i][j]=A[bx-1 +i][j];   
	 A[nx+bx-1+i][j]=A[nx+bx-i][j];
     }
	 
return;
}

/*--------------------------------------------------------------------------*/

#endif
