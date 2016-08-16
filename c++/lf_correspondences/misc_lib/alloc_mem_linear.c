/*****************************************************************************/
/*                                                                           */
/*               Copyright 08/2006 by Dr. Andres Bruhn and,                  */
/*     Faculty of Mathematics and Computer Science, Saarland University,     */
/*                           Saarbruecken, Germany.                          */
/*                                                                           */
/*****************************************************************************/


#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef ALLOC_MEM_LINEAR
#define ALLOC_MEM_LINEAR


/*--------------------------------------------------------------------------*/


int alloc_mem_linear_core

(                   /********************************************************/
 void  **ptr,       /* void pointer on n-dim array - use &((void*)(pointer) */
 long  elem_size,   /* size of each element        - use sizeof(elem)       */
 long  dim,         /* dimension of n-dim array    - use n                  */
 long *dimsize      /* array of dimension sizes    - see example            */
                    /* -----------------------------------------------------*/
                    /* the following example shall help to explain the      */
                    /* source code :   for a[2][3][4][5] ->  dimsize[0]=2   */
                    /*                                       dimsize[1]=3   */
                    /*                                       dimsize[2]=4   */
                    /*                                       dimsize[3]=5   */
)                   /********************************************************/

/* linear allocation of a n-dimensional array  */

{                   /********************************************************/
  long i;           /* counter variable                                     */
  long i_dim;       /* counter variable for dimensions                      */
  long total;       /* total number of elements (finest hierarchy level)    */
  long sum;         /* number of pointer on each hierarchy level            */
  char *totalspace; /* pointer on finest hierarchy level                    */
  char ***ptr_array;/* pointer array on all required pointers               */
                    /********************************************************/

  /*---------------------------- ELEMENT MEMORY ----------------------------*/
//printf(" aloocate core 1 %d \n",dim);
  /**************************************************************************/
  /* counts total number of elements                                        */
  /* -> in our example total = 2*3*4*5 = 120                                */
  /**************************************************************************/
  total = 1;
  for (i_dim=0;i_dim<dim;i_dim++) total*=dimsize[i_dim];
//printf(" aloocate core 2 %d \n",dim);
  /**************************************************************************/
  /* allocates space for all elements                                       */
  /* -> in our example space for 120 elements of elem_size is allocated     */
  /**************************************************************************/
  //printf(" aloocate core 2 %x  %d %d\n",totalspace, total, elem_size);

  //totalspace =  (char*) malloc (total * elem_size);
  totalspace = new char[total * elem_size];
//printf(" aloocate core 3 %x  %d %d\n",totalspace, total, elem_size);
  /**************************************************************************/
  /* if allocation fails                                                    */
  /**************************************************************************/
  if (totalspace==NULL)
  {
      printf(" \n MEMORY ALLOCATION FAILURE \n");
      exit(1);
  }

  /**************************************************************************/
  /* if we have the one dimensional case                                    */
  /**************************************************************************/
  if(dim==1)
  {
      /* set master pointer to allocated space */
      *ptr=((void*)totalspace);
      return(0);
  }


  /**************************************************************************/
  /* if we have more than one dimension                                     */
  /* -> this is the case in our example                                     */
  /**************************************************************************/

  /*---------------------------- POINTER MEMORY ----------------------------*/

  /**************************************************************************/
  /* allocates a pointer for each hierarchy level but the finest.           */
  /* therefore an array of (dim-1) pointers is allocated.                   */
  /* on the finest level we need space for the elements not for pointers    */
  /* this space has already been allocated (totalspace)                     */
  /* -> in our example we want to have:                                     */
  /*       granularity level  expression     type      our memory pointer   */
  /*       coarsest    0      a[2]           pointer   ptr_array[0]         */
  /*                   1      a[2][3]        pointer   ptr_array[1]         */
  /*                   2      a[2][3][4]     pointer   ptr_array[2]         */
  /*       finest      3      a[2][3][4][5]  element   totalspace           */
  /*    therefore we need three pointers on pointer vectors:                */
  /*                                        one for level 0 -> a[2]         */
  /*                                        one for level 1 -> a[2][3]      */
  /*                                        one for level 2 -> a[2][3][4]   */
  /*    ptr_array[0] points on the pointer (vector) on level 0              */
  /*    ptr_array[1] points on the pointer (vector) on level 1              */
  /*    ptr_array[2] points on the pointer (vector) on level 2              */
  /**************************************************************************/
//printf(" aloocate core 4 %d \n",dim);
  //ptr_array = (char***) malloc ((dim-1) * sizeof(char**));
  ptr_array = new char**[(dim-1) * sizeof(char**)];
  /**************************************************************************/
  /* if allocation fails                                                    */
  /**************************************************************************/
  if (ptr_array==NULL)
  {
      printf(" \n MEMORY ALLOCATION FAILURE \n");
      exit(1);
  }
  /**************************************************************************/
  /* allocates space for the pointer vector on each hierarchy level.        */
  /* sum is the number of pointers of this vector on each hierarchy level.  */
  /* sum is computed from coarse to fine (i_dim++)                          */
  /* -> in our example the number of pointers are                           */
  /*    for level 0 :  sum = 1 * 2     =  2                                 */
  /*    for level 1 :  sum = 1 * 2*3   =  6                                 */
  /*    for level 2 :  sum = 1 * 2*3*4 = 24                                 */
  /**************************************************************************/
  sum = 1;
  for (i_dim=0;i_dim<dim-1;i_dim++)
  {
      sum*=dimsize[i_dim];
      //ptr_array[i_dim]= (char**) malloc (sum * sizeof(char*));
      ptr_array[i_dim] = new char*[sum * sizeof(char*)];
      /* if allocation fails */
      if (ptr_array[i_dim]==NULL)
      {
	  printf(" \n MEMORY ALLOCATION FAILURE \n");
	  exit(1);
      }
  }

  /**************************************************************************/
  /* sum is again the number of pointers on each hierarchy level.           */
  /* this time we divide because we count from fine to coarse.              */
  /* -> in our example sum starts with :                                    */
  /*    for level 2   :  sum = 120/5 = 24                                   */
  /**************************************************************************/
  sum=total;
  sum=sum/dimsize[dim-1];

  /**************************************************************************/
  /* sets pointers from the second finest level to the beginning of the     */
  /* element rows on the finest level.                                      */
  /* since pointers point on elements this is a special case!               */
  /* in all other cases pointers point on other pointers.                   */
  /* -> in our example the second finest level is ptr_array[2] (level 2).   */
  /*    the finest level is totalspace (level 3).                           */
  /*    all pointers (24) of ptr_array[2] point on every 5th element with   */
  /*    size elem_size lying in the memory allocated in totalspace.         */
  /*    after this step (ptr_array[2][0..24])[0..4] is accessible           */
  /*    - each 5th element  ->               [0..4]                         */
  /**************************************************************************/
  for(i=0;i<sum;i++)
  {
      ptr_array[dim-2][i]=(char*)(totalspace+(i*dimsize[dim-1])*elem_size);
  }

  /**************************************************************************/
  /* set pointers from level (i_dim) to the beginning of the pointer rows   */
  /* on the next finer level (i_dim+1).                                     */
  /* this is done frome fine to coarse (i_dim--)                            */
  /* -> in our example the third finest level is ptr_array[1] (level 1).    */
  /*    the second finest level is ptr_array[2] (level 2).                  */
  /*    all pointers (6) of ptr_array[1] point on every 4th pointer lying   */
  /*    in the memory allocated in ptr_array[2].                            */
  /*    after this step (ptr_array[1][0..6])[0..3][0..4] is accessible      */
  /*    - each 4th pointer  ->              [0..3]                          */
  /* -> in our example the fourth finest level is ptr_array[0] (level 0).   */
  /*    the third finest level is ptr_array[1] (level 1).                   */
  /*    all pointers (2) of ptr_array[0] point on every 3th pointer lying   */
  /*    in the memory allocated in ptr_array[1].                            */
  /*    after this step (ptr_array[0][0..1])[0..2][0..3][0..4] is accessible*/
  /*    - each 3rd pointer  ->              [0..2]                          */
  /* -> since ptr_array[0] has always the right dimension this can be       */
  /*    reformulated as (ptr_array[0])[0..1][0..2][0..3][0..4].             */
  /*    the pointer hierarchy is complete and ptr_array[0] must become the  */
  /*    master pointer.                                                     */
  /**************************************************************************/


  for(i_dim=dim-3;i_dim>=0;i_dim--)
  {
      sum=sum/dimsize[i_dim+1];
      for(i=0;i<sum;i++)
       ptr_array[i_dim][i]=(char*)(&(ptr_array[i_dim+1][i*dimsize[i_dim+1]]));
  }
  //printf(" aloocate core 2 %d \n",dim);

  /**************************************************************************/
  /* set master pointer to pointer on coarsest level                        */
  /**************************************************************************/
  *ptr=((void*)ptr_array[0]);

  /**************************************************************************/
  /* free ptr_array array                                                   */
  /**************************************************************************/
  free(ptr_array);

  return(0);
}

/*--------------------------------------------------------------------------*/


int alloc_mem_linear

(                   /********************************************************/
 void  **ptr,       /* void pointer on n-dim array - use &((void*)(pointer) */
 long  elem_size,   /* size of each element        - use sizeof(elem)       */
 long  dim,         /* dimension of n-dim array    - use n                  */
 long  args,        /* size of each dimension -> n arguments                */
 ...                /*                                                      */
)                   /********************************************************/

/* user friendly interface for (allows size parameters for each dimension) */

{                   /********************************************************/
  long i;           /* counter variable                                     */
  long *dimsize;    /* array of dimensions                                  */
  va_list ap;       /* object for traversing the arguments                  */
                    /********************************************************/

  /*----------------------------- ARGUMENT LIST ----------------------------*/

  /**************************************************************************/
  /* if number of dimensions smaller or equal than 0                        */
  /**************************************************************************/
  if (dim<=0)
  {
      printf("\n WRONG ARRAY DIMENSIONS \n");
      exit(1);
  }

  /**************************************************************************/
  /* allocate array for dimension sizes                                     */
  /**************************************************************************/
  //dimsize =  (long*) malloc (dim * sizeof(long));
  dimsize = new long[dim *sizeof(long)];

  /**************************************************************************/
  /* if allocation fails                                                    */
  /**************************************************************************/
  if (dimsize==NULL)
  {
      printf(" \n MEMORY ALLOCATION FAILURE \n");
      exit(1);
  }

  /**************************************************************************/
  /* set iterator to begin of argument list                                 */
  /**************************************************************************/
  va_start(ap, args);

  /**************************************************************************/
  /* traversing the argument list                                           */
  /**************************************************************************/
  for(i=0;i<dim;i++)
  {
      /* read argument */
      dimsize[i] = args;
      /* if size of dimension smaller or equal than 0 */
      if (dimsize[i]<=0)
      {
	  printf("\n WRONG ARRAY DIMENSION SIZES \n");
	  exit(1);
      }
      /* next argument */
      args=va_arg(ap, long);
  }

  /**************************************************************************/
  /* close argument list                                                    */
  /**************************************************************************/
  va_end(ap);
//printf(" aloocate inside 1 %d %d\n",elem_size,dim);
  /**************************************************************************/
  /* call core routine with dimsize array                                   */
  /**************************************************************************/
  alloc_mem_linear_core(ptr,elem_size,dim,dimsize);
  //printf(" aloocate inside 2 \n");

  /**************************************************************************/
  /* free dimsize array                                                     */
  /**************************************************************************/
  free(dimsize);

  return(0);
}


/*--------------------------------------------------------------------------*/


int free_mem_linear_core

(                   /********************************************************/
 void  *ptr ,       /* void n-dim array            - use ((void*)(pointer)  */
 long  elem_size,   /* size of each element        - use sizeof(elem)       */
 long  dim,         /* dimension of n-dim array    - use n                  */
 long  *dimsize     /* array of dimension sizes    - see example            */
                    /* -----------------------------------------------------*/
                    /* the following example shall help to explain the      */
                    /* source code :   for a[2][3][4] ->  dimsize[0]=2      */
                    /*                                    dimsize[1]=3      */
                    /*                                    dimsize[2]=4      */
)                   /********************************************************/

/* free n-dimensional array allocated by alloc_mem_linear_core  */

{                   /********************************************************/
  long i_dim;       /* counter variable for dimensions                      */
  long j_dim;       /* counter variable for dimensions                      */
  char** tmpptr1;   /* pointer to recursively free all allocated structures */
                    /********************************************************/


  /**************************************************************************/
  /* free all pointer levels from fine to coarse (j_dim--).                 */
  /* j_dim is the number of required recursions in order to obtain the      */
  /* pointer on the desired level.                                          */
  /* -> in our example the first pointer that is freed is a[0][0][0].       */
  /*    this pointer points on the first element of our n-dim array and     */
  /*    frees all the memory used by the elements.                          */
  /*    thereafter all pointer that are freed (a[0][0] and a[0) belong to   */
  /*    the pointer hierarchy.                                              */
  /**************************************************************************/
  for(j_dim=dim-1;j_dim>=1;j_dim--)
  {
      /* set pointer to coarsest level */
      tmpptr1=(char**)(ptr);
      /* performing j_dim recursions to obtain the */
      /* pointer on the desired level              */
      for(i_dim=0;i_dim<j_dim;i_dim++)
      {
	  /* recursion step */
	  tmpptr1=(char**)(tmpptr1[0]);
      }
      /* free pointer */
      free((char*)tmpptr1);
  }


  /**************************************************************************/
  /* free master pointer on coarsest hierarchy level                        */
  /* -> finally the master pointer is freed                                 */
  /**************************************************************************/
  free(ptr);

  return(0);
}

/*--------------------------------------------------------------------------*/


int free_mem_linear

(                   /********************************************************/
 void  *ptr ,       /* void n-dim array            - use ((void*)(pointer)  */
 long  elem_size,   /* size of each element        - use sizeof(elem)       */
 long  dim,         /* dimension of n-dim array    - use n                  */
 long  args,        /* size of each dimension -> n arguments                */
 ...                /*                                                      */
)                   /********************************************************/

/* user friendly interface for (allows size parameters for each dimension) */

{                   /********************************************************/
  long i;           /* counter variable                                     */
  long *dimsize;    /* array of dimensions                                  */
  va_list ap;       /* object for traversing the arguments                  */
                    /********************************************************/


  /**************************************************************************/
  /* if number of dimensions smaller or equal than 0                        */
  /**************************************************************************/
  if (dim<=0)
  {
      printf("\n WRONG ARRAY DIMENSIONS \n");
      exit(1);
  }

  /**************************************************************************/
  /* allocate array for dimension sizes                                     */
  /**************************************************************************/
  //dimsize =  (long*) malloc (dim * sizeof(long));
  dimsize = new long[dim * sizeof(long)];
  /**************************************************************************/
  /* if allocation fails                                                    */
  /**************************************************************************/
  if (dimsize==NULL)
  {
      printf(" \n MEMORY ALLOCATION FAILURE \n");
      exit(1);
  }

  /**************************************************************************/
  /* set iterator to begin of argument list                                 */
  /**************************************************************************/
  va_start(ap, args);

  /**************************************************************************/
  /* traversing the argument list                                           */
  /**************************************************************************/
  for(i=0;i<dim;i++)
  {
      /* read argument */
      dimsize[i] = args;
      /* if size of dimension smaller or equal than 0 */
      if (dimsize[i]<=0)
      {
	  printf("\n WRONG ARRAY DIMENSION SIZES \n");
      }
      /* next argument */
      args=va_arg(ap, long);
  }

  /**************************************************************************/
  /* close argument list                                                    */
  /**************************************************************************/
  va_end(ap);

  /**************************************************************************/
  /* call core routine with dimsize array                                   */
  /**************************************************************************/
  free_mem_linear_core(ptr,elem_size,dim,dimsize);

  /**************************************************************************/
  /* free dimsize array                                                     */
  /**************************************************************************/
  free(dimsize);

  return(0);
}


/*--------------------------------------------------------------------------*/

#endif
