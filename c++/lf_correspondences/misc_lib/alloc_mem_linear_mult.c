/*****************************************************************************/
/*                                                                           */
/*                   Copyright 08/2006 by Dr. Andres Bruhn,                  */
/*     Faculty of Mathematics and Computer Science, Saarland University,     */
/*                           Saarbruecken, Germany.                          */
/*                                                                           */
/*****************************************************************************/


#ifndef ALLOC_MEM_LINEAR_MULT
#define ALLOC_MEM_LINEAR_MULT

#include <stdarg.h>
#include "alloc_mem_linear.c"


/*--------------------------------------------------------------------------*/

int ALLOC_VECTOR

(                   /********************************************************/
 int   n,           /* number of vectors                                    */
 int   nx,          /* size in x direction                                  */
 float **args,      /* pointer for each vector -> n arguments               */
 ...                /*                                                      */
)                   /********************************************************/

{                      /*****************************************************/
  float ***ptr_array;  /* auxiliary pointer for arguments                   */
  float **call_ptr;    /* auxiliary pointer for function call               */
  int   i;             /* counter variable                                  */
  va_list ap;          /* object for traversing the arguments               */
                       /*****************************************************/

  /*----------------------------- ARGUMENT LIST ----------------------------*/

  /**************************************************************************/
  /* if number of vectors smaller or equal than 0                           */
  /**************************************************************************/
  if (n<=0)
  {
      printf("\n ALLOC_VECTOR : WRONG NUMBER OF VECTORS \n");
      exit(1);
  }

  /**************************************************************************/
  /* allocate array for arguments                                           */
  /**************************************************************************/
  //ptr_array =  (float***) malloc (n * sizeof(float**));
  ptr_array = new float**[n * sizeof(float**)];
  /**************************************************************************/
  /* if allocation fails                                                    */
  /**************************************************************************/
  if (ptr_array==NULL)
  {
      printf(" \n ALLOC_VECTOR : MEMORY ALLOCATION FAILURE \n");
      exit(1);
  }

  /**************************************************************************/
  /* set iterator to begin of argument list                                 */
  /**************************************************************************/
  va_start(ap, args);

  /**************************************************************************/
  /* traversing the argument list                                           */
  /**************************************************************************/
  for(i=0;i<n;i++)
  {
      /* read argument */
      ptr_array[i] = args;
      /* if NULL POINTER argument */
      if (ptr_array[i]==NULL)
      {
	  printf("\n ALLOC_VECTOR : NULL POINTER ADDRESS EXCEPTION\n");
	  exit(1);
      }
      /* next argument */
      args=va_arg(ap, float**);
  }

  /**************************************************************************/
  /* close argument list                                                    */
  /**************************************************************************/
  va_end(ap);

  /**************************************************************************/
  /* call core routine with call_ptr                                        */
  /**************************************************************************/
  alloc_mem_linear((void**)(&call_ptr),sizeof(float), 2, n, nx);

  /**************************************************************************/
  /* redistribute allocated memory                                          */
  /**************************************************************************/
  for(i=0;i<n;i++)
  {
      *(ptr_array[i])=call_ptr[i];
  }

  /**************************************************************************/
  /* free memory                                                            */
  /**************************************************************************/
  free(ptr_array);
  free(call_ptr);

  return(0);
}


/*--------------------------------------------------------------------------*/


int ALLOC_MATRIX

(                   /********************************************************/
 int   n,           /* number of matrices                                   */
 int   nx,          /* size in x direction                                  */
 int   ny,          /* size in y direction                                  */
 float ***args,     /* pointer for each matrix -> n arguments               */
 ...                /*                                                      */
)                   /********************************************************/

{                      /*****************************************************/
  float ****ptr_array; /* auxiliary pointer for arguments                   */
  float ***call_ptr;   /* auxiliary pointer for function call               */
  int   i;             /* counter variable                                  */
  va_list ap;          /* object for traversing the arguments               */
                       /*****************************************************/

  /*----------------------------- ARGUMENT LIST ----------------------------*/

  /**************************************************************************/
  /* if number of matrices smaller or equal than 0                          */
  /**************************************************************************/
  if (n<=0)
  {
      printf("\n ALLOC_MATRIX : WRONG NUMBER OF MATRICES \n");
      exit(1);
  }

  /**************************************************************************/
  /* allocate array for arguments                                           */
  /**************************************************************************/
  //ptr_array =  (float****) malloc (n * sizeof(float***));
  ptr_array = new float***[n * sizeof(float***)];
  /**************************************************************************/
  /* if allocation fails                                                    */
  /**************************************************************************/
  if (ptr_array==NULL)
  {
      printf(" \n ALLOC_MATRIX : MEMORY ALLOCATION FAILURE \n");
      exit(1);
  }

  /**************************************************************************/
  /* set iterator to begin of argument list                                 */
  /**************************************************************************/
  va_start(ap, args);

  /**************************************************************************/
  /* traversing the argument list                                           */
  /**************************************************************************/
  for(i=0;i<n;i++)
  {
      /* read argument */
      ptr_array[i] = args;
      /* if NULL POINTER argument */
      if (ptr_array[i]==NULL)
      {
	  printf("\n ALLOC_MATRIX : NULL POINTER ADDRESS EXCEPTION\n");
	  exit(1);
      }
      /* next argument */
      args=va_arg(ap, float***);
  }

  /**************************************************************************/
  /* close argument list                                                    */
  /**************************************************************************/
  va_end(ap);

  /**************************************************************************/
  /* call core routine with call_ptr                                        */
  /**************************************************************************/
  alloc_mem_linear((void**)(&call_ptr),sizeof(float), 3, n, nx, ny);

  /**************************************************************************/
  /* redistribute allocated memory                                          */
  /**************************************************************************/

  for(i=0;i<n;i++)
  {
      *(ptr_array[i])=call_ptr[i];
  }

  /**************************************************************************/
  /* free memory                                                            */
  /**************************************************************************/
  free(ptr_array);
  free(call_ptr);

  return(0);
}

/*--------------------------------------------------------------------------*/


int ALLOC_CUBIX

(                   /********************************************************/
 int   n,           /* number of cubices                                    */
 int   nx,          /* size in x direction                                  */
 int   ny,          /* size in y direction                                  */
 int   nz,          /* size in z direction                                  */
 float ****args,    /* pointer for each cubix -> n arguments                */
 ...                /*                                                      */
)                   /********************************************************/

{                      /*****************************************************/
  float *****ptr_array;/* auxiliary pointer for arguments                   */
  float ****call_ptr;  /* auxiliary pointer for function call               */
  int   i;             /* counter variable                                  */
  va_list ap;          /* object for traversing the arguments               */
                       /*****************************************************/

  /*----------------------------- ARGUMENT LIST ----------------------------*/

  /**************************************************************************/
  /* if number of cubiices smaller or equal than 0                          */
  /**************************************************************************/
  if (n<=0)
  {
      printf("\n ALLOC_CUBIX : WRONG NUMBER OF CUBICES \n");
      exit(1);
  }

  //printf(" aloocate step 0 \n");
  /**************************************************************************/
  /* allocate array for arguments                                           */
  /**************************************************************************/
  //ptr_array =  (float*****) malloc (n * sizeof(float****));
  ptr_array = new float****[n * sizeof(float****)];
  /**************************************************************************/
  /* if allocation fails                                                    */
  /**************************************************************************/
  if (ptr_array==NULL)
  {
      printf(" \n ALLOC_CUBIX : MEMORY ALLOCATION FAILURE \n");
      exit(1);
  }

  //printf(" aloocate step 1 \n");
  /**************************************************************************/
  /* set iterator to begin of argument list                                 */
  /**************************************************************************/
  va_start(ap, args);

  /**************************************************************************/
  /* traversing the argument list                                           */
  /**************************************************************************/
  for(i=0;i<n;i++)
  {
      /* read argument */
      ptr_array[i] = args;
      /* if NULL POINTER argument */
      if (ptr_array[i]==NULL)
      {
	  printf("\n ALLOC_CUBIX : NULL POINTER ADDRESS EXCEPTION\n");
	  exit(1);
      }
      /* next argument */
      args=va_arg(ap, float****);
  }

  /**************************************************************************/
  /* close argument list                                                    */
  /**************************************************************************/
  va_end(ap);
  //printf(" aloocate step 2 \n");
  /**************************************************************************/
  /* call core routine with call_ptr                                        */
  /**************************************************************************/
  alloc_mem_linear((void**)(&call_ptr),sizeof(float), 4, n, nx, ny, nz);
//printf(" aloocate step 3 \n");
  /**************************************************************************/
  /* redistribute allocated memory                                          */
  /**************************************************************************/

  for(i=0;i<n;i++)
  {
      *(ptr_array[i])=call_ptr[i];
  }

  /**************************************************************************/
  /* free memory                                                            */
  /**************************************************************************/
  free(ptr_array);
  free(call_ptr);

  return(0);
}

/*--------------------------------------------------------------------------*/


int ALLOC_QUADRIX

(                   /********************************************************/
 int   n,           /* number of quadrices                                  */
 int   na,          /* size in a direction                                  */
 int   nx,          /* size in x direction                                  */
 int   ny,          /* size in y direction                                  */
 int   nz,          /* size in z direction                                  */
 float *****args,   /* pointer for each quadrix -> n arguments              */
 ...                /*                                                      */
)                   /********************************************************/

{                       /****************************************************/
  float ******ptr_array;/* auxiliary pointer for arguments                  */
  float *****call_ptr;  /* auxiliary pointer for function call              */
  int   i;              /* counter variable                                 */
  va_list ap;           /* object for traversing the arguments              */
                        /****************************************************/

  /*----------------------------- ARGUMENT LIST ----------------------------*/

  /**************************************************************************/
  /* if number of quadrices smaller or equal than 0                         */
  /**************************************************************************/
  if (n<=0)
  {
      printf("\n ALLOC_QUADRIX : WRONG NUMBER OF QUADRICES \n");
      exit(1);
  }

  /**************************************************************************/
  /* allocate array for arguments                                           */
  /**************************************************************************/
  //ptr_array =  (float******) malloc (n * sizeof(float*****));
  ptr_array =  new float*****[n * sizeof(float*****)];
  /**************************************************************************/
  /* if allocation fails                                                    */
  /**************************************************************************/
  if (ptr_array==NULL)
  {
      printf(" \n ALLOC_QUADRIX : MEMORY ALLOCATION FAILURE \n");
      exit(1);
  }

  /**************************************************************************/
  /* set iterator to begin of argument list                                 */
  /**************************************************************************/
  va_start(ap, args);

  /**************************************************************************/
  /* traversing the argument list                                           */
  /**************************************************************************/
  for(i=0;i<n;i++)
  {
      /* read argument */
      ptr_array[i] = args;
      /* if NULL POINTER argument */
      if (ptr_array[i]==NULL)
      {
	  printf("\n ALLOC_QUADRIX : NULL POINTER ADDRESS EXCEPTION\n");
	  exit(1);
      }
      /* next argument */
      args=va_arg(ap, float*****);
  }

  /**************************************************************************/
  /* close argument list                                                    */
  /**************************************************************************/
  va_end(ap);

  /**************************************************************************/
  /* call core routine with call_ptr                                        */
  /**************************************************************************/
  alloc_mem_linear((void**)(&call_ptr),sizeof(float), 5, n, na, nx, ny, nz);

  /**************************************************************************/
  /* redistribute allocated memory                                          */
  /**************************************************************************/

  for(i=0;i<n;i++)
  {
      *(ptr_array[i])=call_ptr[i];
  }

  /**************************************************************************/
  /* free memory                                                            */
  /**************************************************************************/
  free(ptr_array);
  free(call_ptr);

  return(0);
}

/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/


int FREE_VECTOR

(                   /********************************************************/
 int   n,           /* number of vectors                                    */
 int   nx,          /* size in x direction                                  */
 float *args,       /* vectors -> n arguments                               */
 ...                /*                                                      */
)                   /********************************************************/

{                      /*****************************************************/
  float **ptr_array;   /* auxiliary pointer for arguments                   */
  int   i;             /* counter variable                                  */
  va_list ap;          /* object for traversing the arguments               */
                       /*****************************************************/

  /*----------------------------- ARGUMENT LIST ----------------------------*/

  /**************************************************************************/
  /* if number of vectors smaller or equal than 0                           */
  /**************************************************************************/
  if (n<=0)
  {
      printf("\n FREE_VECTOR : WRONG NUMBER OF VECTORS \n");
      exit(1);
  }

  /**************************************************************************/
  /* allocate array for arguments                                           */
  /**************************************************************************/
  //ptr_array =  (float**) malloc (n * sizeof(float*));
  ptr_array =  new float*[n * sizeof(float*)];

  /**************************************************************************/
  /* if allocation fails                                                    */
  /**************************************************************************/
  if (ptr_array==NULL)
  {
      printf(" \n FREE_VECTOR : MEMORY ALLOCATION FAILURE \n");
      exit(1);
  }

  /**************************************************************************/
  /* set iterator to begin of argument list                                 */
  /**************************************************************************/
  va_start(ap, args);

  /**************************************************************************/
  /* traversing the argument list                                           */
  /**************************************************************************/
  for(i=0;i<n;i++)
  {
      /* read argument */
      ptr_array[i] = args;
      /* if NULL POINTER argument */
      if (ptr_array[i]==NULL)
      {
	  printf("\n FREE_VECTOR : NULL POINTER ADDRESS EXCEPTION\n");
	  exit(1);
      }
      /* next argument */
      args=va_arg(ap, float*);
  }

  /**************************************************************************/
  /* close argument list                                                    */
  /**************************************************************************/
  va_end(ap);


  /**************************************************************************/
  /* free memory of main vector                                             */
  /**************************************************************************/
  free_mem_linear((void*)(ptr_array[0]),sizeof(float), 1, nx);

  /**************************************************************************/
  /* free memory                                                            */
  /**************************************************************************/
  free(ptr_array);

  return(0);
}

/*--------------------------------------------------------------------------*/


int FREE_MATRIX

(                   /********************************************************/
 int   n,           /* number of matrices                                   */
 int   nx,          /* size in x direction                                  */
 int   ny,          /* size in y direction                                  */
 float **args,      /* matrices -> n arguments                              */
 ...                /*                                                      */
)                   /********************************************************/

{                      /*****************************************************/
  float ***ptr_array;  /* auxiliary pointer for arguments                   */
  int   i;             /* counter variable                                  */
  va_list ap;          /* object for traversing the arguments               */
                       /*****************************************************/

  /*----------------------------- ARGUMENT LIST ----------------------------*/

  /**************************************************************************/
  /* if number of matrices is smaller or equal than 0                       */
  /**************************************************************************/
  if (n<=0)
  {
      printf("\n FREE_MATRIX : WRONG NUMBER OF MATRICES \n");
      exit(1);
  }

  /**************************************************************************/
  /* allocate array for arguments                                           */
  /**************************************************************************/
  //ptr_array =  (float***) malloc (n * sizeof(float**));
  ptr_array =  new float**[(n * sizeof(float**))];

  /**************************************************************************/
  /* if allocation fails                                                    */
  /**************************************************************************/
  if (ptr_array==NULL)
  {
      printf(" \n FREE_MATRIX : MEMORY ALLOCATION FAILURE \n");
      exit(1);
  }

  /**************************************************************************/
  /* set iterator to begin of argument list                                 */
  /**************************************************************************/
  va_start(ap, args);

  /**************************************************************************/
  /* traversing the argument list                                           */
  /**************************************************************************/
  for(i=0;i<n;i++)
  {
      /* read argument */
      ptr_array[i] = args;
      /* if NULL POINTER argument */
      if (ptr_array[i]==NULL)
      {
	  printf("\n FREE_MATRIX : NULL POINTER ADDRESS EXCEPTION\n");
	  exit(1);
      }
      /* next argument */
      args=va_arg(ap, float**);
  }

  /**************************************************************************/
  /* close argument list                                                    */
  /**************************************************************************/
  va_end(ap);


  /**************************************************************************/
  /* free memory of main matrix                                             */
  /**************************************************************************/
  free_mem_linear((void*)(ptr_array[0]),sizeof(float), 2, nx, ny);

  /**************************************************************************/
  /* free memory                                                            */
  /**************************************************************************/
  free(ptr_array);

  return(0);
}

/*--------------------------------------------------------------------------*/


int FREE_CUBIX

(                   /********************************************************/
 int   n,           /* number of cubices                                    */
 int   nx,          /* size in x direction                                  */
 int   ny,          /* size in y direction                                  */
 int   nz,          /* size in z direction                                  */
 float ***args,     /* cubices -> n arguments                               */
 ...                /*                                                      */
)                   /********************************************************/

{                      /*****************************************************/
  float ****ptr_array; /* auxiliary pointer for arguments                   */
  int   i;             /* counter variable                                  */
  va_list ap;          /* object for traversing the arguments               */
                       /*****************************************************/

  /*----------------------------- ARGUMENT LIST ----------------------------*/

  /**************************************************************************/
  /* if number of cubices is smaller or equal than 0                        */
  /**************************************************************************/
  if (n<=0)
  {
      printf("\n FREE_CUBIX : WRONG NUMBER OF CUBICES \n");
      exit(1);
  }

  /**************************************************************************/
  /* allocate array for arguments                                           */
  /**************************************************************************/
  //ptr_array =  (float****) malloc (n * sizeof(float***));
  ptr_array =  new float***[(n * sizeof(float***))];

  /**************************************************************************/
  /* if allocation fails                                                    */
  /**************************************************************************/
  if (ptr_array==NULL)
  {
      printf(" \n FREE_CUBIX : MEMORY ALLOCATION FAILURE \n");
      exit(1);
  }

  /**************************************************************************/
  /* set iterator to begin of argument list                                 */
  /**************************************************************************/
  va_start(ap, args);

  /**************************************************************************/
  /* traversing the argument list                                           */
  /**************************************************************************/
  for(i=0;i<n;i++)
  {
      /* read argument */
      ptr_array[i] = args;
      /* if NULL POINTER argument */
      if (ptr_array[i]==NULL)
      {
	  printf("\n FREE_CUBIX : NULL POINTER ADDRESS EXCEPTION\n");
	  exit(1);
      }
      /* next argument */
      args=va_arg(ap, float***);
  }

  /**************************************************************************/
  /* close argument list                                                    */
  /**************************************************************************/
  va_end(ap);


  /**************************************************************************/
  /* free memory of main cubix                                              */
  /**************************************************************************/
  free_mem_linear((void*)(ptr_array[0]),sizeof(float), 3, nx, ny, nz);

  /**************************************************************************/
  /* free memory                                                            */
  /**************************************************************************/
  free(ptr_array);

  return(0);
}

/*--------------------------------------------------------------------------*/


int FREE_QUADRIX

(                   /********************************************************/
 int   n,           /* number of quadrices                                  */
 int   na,          /* size in a direction                                  */
 int   nx,          /* size in x direction                                  */
 int   ny,          /* size in y direction                                  */
 int   nz,          /* size in z direction                                  */
 float ****args,    /* quadrices -> n arguments                             */
 ...                /*                                                      */
)                   /********************************************************/

{                       /****************************************************/
  float *****ptr_array; /* auxiliary pointer for arguments                  */
  int   i;              /* counter variable                                 */
  va_list ap;           /* object for traversing the arguments              */
                        /****************************************************/

  /*----------------------------- ARGUMENT LIST ----------------------------*/

  /**************************************************************************/
  /* if number of quadrices is smaller or equal than 0                      */
  /**************************************************************************/
  if (n<=0)
  {
      printf("\n FREE_QUADRIX : WRONG NUMBER OF QUADRICES \n");
      exit(1);
  }

  /**************************************************************************/
  /* allocate array for arguments                                           */
  /**************************************************************************/
  //ptr_array =  (float*****) malloc (n * sizeof(float****));
  ptr_array =  new float****[(n * sizeof(float****))];

  /**************************************************************************/
  /* if allocation fails                                                    */
  /**************************************************************************/
  if (ptr_array==NULL)
  {
      printf(" \n FREE_QUADRIX : MEMORY ALLOCATION FAILURE \n");
      exit(1);
  }

  /**************************************************************************/
  /* set iterator to begin of argument list                                 */
  /**************************************************************************/
  va_start(ap, args);

  /**************************************************************************/
  /* traversing the argument list                                           */
  /**************************************************************************/
  for(i=0;i<n;i++)
  {
      /* read argument */
      ptr_array[i] = args;
      /* if NULL POINTER argument */
      if (ptr_array[i]==NULL)
      {
	  printf("\n FREE_QUADRIX : NULL POINTER ADDRESS EXCEPTION\n");
	  exit(1);
      }
      /* next argument */
      args=va_arg(ap, float****);
  }

  /**************************************************************************/
  /* close argument list                                                    */
  /**************************************************************************/
  va_end(ap);


  /**************************************************************************/
  /* free memory of main quadrix                                            */
  /**************************************************************************/
  free_mem_linear((void*)(ptr_array[0]),sizeof(float), 4, na, nx, ny, nz);

  /**************************************************************************/
  /* free memory                                                            */
  /**************************************************************************/
  free(ptr_array);

  return(0);
}


/*--------------------------------------------------------------------------*/

#endif
