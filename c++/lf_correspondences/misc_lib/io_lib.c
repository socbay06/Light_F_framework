#ifndef IO_LIB_INCLUDED
#define IO_LIB_INCLUDED

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc_mem_linear.c"
#include "color_lib.c"
#include "matrix_lib.c"
#include "funct_lib.c"



/*--------------------------------------------------------------------------*/

void read_pgm_header
(
                            /************************************************/
    char  *filename,        /* in     : name of pgm file                    */
    long  *header_end_pos,  /* in+out : position of header end              */
    int   *nx,              /* out    : size in x-direction                 */
    int   *ny               /* out    : size in y-direction                 */
                            /************************************************/
)
/* reads PGM header */
{
                  /**********************************************************/
    char row[80]; /* read buffer                                            */
    FILE *file;   /* file pointer                                           */
                  /**********************************************************/

    printf("\n Trying to read file %s ...",filename);

    /* try to open file */
    file = fopen(filename,"r");

    /* if file not found */
    if (file==NULL)
    {
	printf("... FAILED");
	printf("\n\n PROGRAM ABORTED !!! \n\n");
	exit(0);
    }

    /* read header */
    fgets (row, 300, file);
    fgets (row, 300, file);
    while (row[0]=='#') fgets(row, 300, file);

    /* read image dimensions */
    sscanf (row, "%d %d", nx, ny);
    fgets (row, 300, file);

    /* determine end of header */
    *header_end_pos=ftell(file);

    /* close file */
    fclose(file);
}

/* ------------------------------------------------------------------------- */


void read_pgm_data
(
                            /************************************************/
    char  *filename,        /* in   : name of PGM file                      */
    long  data_start_pos,   /* in   : position of data start                */
    float **u,              /* out  : image                                 */
    int   nx,               /* in   : size in x-direction                   */
    int   ny,               /* in   : size in y-direction                   */
    int   bx,               /* in   : boundary in x-direction               */
    int   by                /* in   : boundary in y-direction               */
                            /************************************************/
)

/* reads PGM data */

{
                  /**********************************************************/
    int i,j;      /* loop variables                                         */
    FILE *file;   /* file pointer                                           */
                  /**********************************************************/

    /* open file */
    file = fopen(filename,"r");

    /* start at begin of data */
    fseek(file,data_start_pos,SEEK_SET);

    /* read image data */
    for (j=by; j<ny+by; j++)
	for (i=bx; i<nx+bx; i++)
	{
	    u[i][j] = (float) getc (file);
	}

    /* close file */
    fclose(file);

    printf("... SUCCESS");
}


/*--------------------------------------------------------------------------*/


void read_ppm_header
(
                            /************************************************/
    char  *filename,        /* in     : name of pgm file                    */
    long  *header_end_pos,  /* in+out : position of header end              */
    int   *nx,              /* out    : size in x-direction                 */
    int   *ny               /* out    : size in y-direction                 */
                            /************************************************/
)
/* reads PGM header */
{
                  /**********************************************************/
    char row[80]; /* read buffer                                            */
    FILE *file;   /* file pointer                                           */
                  /**********************************************************/

    printf("\n Trying to read file %s ...",filename);

    /* try to open file */
    file = fopen(filename,"r");

    /* if file not found */
    if (file==NULL)
    {
	printf("... FAILED");
	printf("\n\n PROGRAM ABORTED !!! \n\n");
	exit(0);
    }

    /* read header */
    fgets (row, 300, file);
    fgets (row, 300, file);
    while (row[0]=='#') fgets(row, 300, file);

    /* read image dimensions */
    sscanf (row, "%d %d", nx, ny);
    fgets (row, 300, file);

    /* determine end of header */
    *header_end_pos=ftell(file);

    /* close file */
    fclose(file);
}

/* ------------------------------------------------------------------------- */


void read_ppm_data
(
                            /************************************************/
    char  *filename,        /* in   : name of PGM file                      */
    long  data_start_pos,   /* in   : position of data start                */
    float ***rgb,           /* out  : RGB image                             */
    int   nx,               /* in   : size in x-direction                   */
    int   ny,               /* in   : size in y-direction                   */
    int   bx,               /* in   : boundary in x-direction               */
    int   by                /* in   : boundary in y-direction               */
                            /************************************************/
)

/* reads PGM data */

{
                  /**********************************************************/
    int  i,j;     /* loop variables                                         */
    FILE *file;   /* file pointer                                           */
                  /**********************************************************/

    /* open file */
    file = fopen(filename,"r");

    /* start at begin of data */
    fseek(file,data_start_pos,SEEK_SET);

    /* read image data */
    for (j=by; j<ny+by; j++)
	for (i=bx; i<nx+bx; i++)
	{
	    rgb[i][j][0] = (float) getc (file);
	    rgb[i][j][1] = (float) getc (file);
	    rgb[i][j][2] = (float) getc (file);
	}

    /* close file */
    fclose(file);

    printf("... SUCCESS");
}

/* ------------------------------------------------------------------------- */


void read_ppm_data_channelwise
(
                            /************************************************/
    char   *filename,       /* in   : name of PGM file                      */
    long   data_start_pos,  /* in   : position of data start                */
    float  **r,             /* out  : R channel of RGB image                */
    float  **g,             /* out  : G channel of RGB image                */
    float  **b,             /* out  : B channel of RGB image                */
    int    nx,              /* in   : size in x-direction                   */
    int    ny,              /* in   : size in y-direction                   */
    int    bx,              /* in   : boundary in x-direction               */
    int    by               /* in   : boundary in y-direction               */
                            /************************************************/
)

/* reads PGM data */

{
                  /**********************************************************/
    int i,j;      /* loop variables                                         */
    FILE *file;   /* file pointer                                           */
                  /**********************************************************/

    /* open file */
    file = fopen(filename,"r");

    /* start at begin of data */
    fseek(file,data_start_pos,SEEK_SET);

    /* read image data */
    for (j=by; j<ny+by; j++)
	for (i=bx; i<nx+bx; i++)
	{
	    r[i][j] = (float) getc (file);
	    g[i][j] = (float) getc (file);
	    b[i][j] = (float) getc (file);
	}

    /* close file */
    fclose(file);

    printf("... SUCCESS");
}

/* ------------------------------------------------------------------------- */

void write_pgm_blank_header
(
                    /********************************************************/
    char *filename, /* in  : file name                                      */
    int  nx,        /* in  : size in x-direction                            */
    int  ny         /* in  : size in y-direction                            */
                    /********************************************************/
)

/* writes blank PGM header */

{

                  /**********************************************************/
    FILE *file;   /* file pointer                                           */
                  /**********************************************************/

    printf("\n Trying to write file %s ...",filename);

    /* try to open file */
    file = fopen (filename, "w");

    /* if file not found -> break */
    if (file==NULL)
    {
	printf("... FAILED");
	printf("\n\n PROGRAM ABORTED !!! \n\n");
	exit(0);
    }

    /* write out header */
    fprintf (file, "P5 \n");
    fprintf (file, "%d %d \n255\n", (int)nx, (int)ny);

    /* close file */
    fclose(file);
}


/* ------------------------------------------------------------------------- */


void write_pgm_data
(
                            /************************************************/
    char  *filename,        /* in : file name                               */
    float **u,              /* in : image                                   */
    int   nx,               /* in : size in x-direction                     */
    int   ny,               /* in : size in y-direction                     */
    int   bx,               /* in : boundary in x-direction                 */
    int   by                /* in : boundary in y-direction                 */
                            /************************************************/
)

/* writes PGM data */

{
                        /****************************************************/
    FILE *file;         /* file pointer                                     */
    int i,j;            /* loop variables                                   */
    float help;         /* tmp variable                                     */
    unsigned char byte; /* variable for conversion                          */
                        /****************************************************/

    /* open file */
    file = fopen (filename, "a");

    /* write image data */
    for (j=by; j<ny+by; j++)
	for (i=bx; i<nx+bx; i++)
	{
	    help = u[i][j];
	    if (help < 0.0)
		byte = (unsigned char)(0.0);
	    else if (help > 255.0)
		byte = (unsigned char)(255.0);
	    else
		byte = (unsigned char)(help);
	    fwrite (&byte, sizeof(unsigned char), 1, file);
	}

    /* close file */
    fclose(file);

    printf("... SUCCESS");
}

/* ------------------------------------------------------------------------- */

void write_ppm_blank_header
(
                            /************************************************/
    char *filename,         /* in  : file name                              */
    int  nx,                /* in  : size in x-direction                    */
    int  ny                 /* in  : size in y-direction                    */
                            /************************************************/
)

/* writes blank PGM header */

{
                  /**********************************************************/
    FILE *file;   /* file pointer                                           */
                  /**********************************************************/

    printf("\n Trying to write file %s ...",filename);

    /* try to open file */
    file = fopen (filename, "w");

    /* if file not found -> break */
    if (file==NULL)
    {
	printf("... FAILED");
	printf("\n\n PROGRAM ABORTED !!! \n\n");
	exit(0);
    }

    /* write out header */
    fprintf (file, "P6 \n");
    fprintf (file, "%d %d \n255\n", (int)nx, (int)ny);

    /* close file */
    fclose(file);
}


/* ------------------------------------------------------------------------- */


void write_ppm_data
(
                            /************************************************/
    char  *filename,        /* in : file name                               */
    float ***rgb,           /* in : RGB image                               */
    int   nx,               /* in : size in x-direction                     */
    int   ny,               /* in : size in y-direction                     */
    int   bx,               /* in : boundary in x-direction                 */
    int   by                /* in : boundary in y-direction                 */
                            /************************************************/
)

/* writes PGM data */

{

                        /*****************************************************/
    FILE *file ;        /* file pointer                                      */
    int i,j;            /* loop variables                                    */
    float help;         /* tmp variable                                      */
    unsigned char byte; /* variable for conversion                           */
                        /*****************************************************/

    /* open file */
    file = fopen (filename, "a");

    /* write image data */
    for (j=by; j<ny+by; j++)
	for (i=bx; i<nx+bx; i++)
	{
	    help = rgb[i][j][0];
	    if (help < 0.0)
		byte = (unsigned char)(0.0);
	    else if (help > 255.0)
		byte = (unsigned char)(255.0);
	    else
		byte = (unsigned char)(help);
	    fwrite (&byte, sizeof(unsigned char), 1, file);

	    help = rgb[i][j][1];
	    if (help < 0.0)
		byte = (unsigned char)(0.0);
	    else if (help > 255.0)
		byte = (unsigned char)(255.0);
	    else
		byte = (unsigned char)(help);
	    fwrite (&byte, sizeof(unsigned char), 1, file);

	    help = rgb[i][j][2];
	    if (help < 0.0)
		byte = (unsigned char)(0.0);
	    else if (help > 255.0)
		byte = (unsigned char)(255.0);
	    else
		byte = (unsigned char)(help);
	    fwrite (&byte, sizeof(unsigned char), 1, file);
	}

    /* close file */
    fclose(file);

    printf("... SUCCESS");
}

/* ------------------------------------------------------------------------- */


void write_ppm_data_channelwise
(
                            /************************************************/
    char   *filename,       /* in   : name of PGM file                      */
    float **r,              /* out  : R channel of RGB image                */
    float **g,              /* out  : G channel of RGB image                */
    float **b,              /* out  : B channel of RGB image                */
    int  nx,                /* in   : size in x-direction                   */
    int  ny,                /* in   : size in y-direction                   */
    int  bx,                /* in   : boundary in x-direction               */
    int  by                 /* in   : boundary in y-direction               */
                            /************************************************/
)

/* reads PPM data channelwise*/

{
                        /****************************************************/
    FILE *file;         /* file pointer                                     */
    int i,j;            /* loop variables                                   */
    float help;         /* tmp variable                                     */
    unsigned char byte; /* variable for conversion                          */
                        /****************************************************/

    /* open file */
    file = fopen (filename, "a");


    /* write image data */
    for (j=by; j<ny+by; j++)
	for (i=bx; i<nx+bx; i++)
	{
	    help = r[i][j];
	    if (help < 0.0)
		byte = (unsigned char)(0.0);
	    else if (help > 255.0)
		byte = (unsigned char)(255.0);
	    else
		byte = (unsigned char)(help);
	    fwrite (&byte, sizeof(unsigned char), 1, file);

	    help = g[i][j];
	    if (help < 0.0)
		byte = (unsigned char)(0.0);
	    else if (help > 255.0)
		byte = (unsigned char)(255.0);
	    else
		byte = (unsigned char)(help);
	    fwrite (&byte, sizeof(unsigned char), 1, file);

	    help = b[i][j];
	    if (help < 0.0)
		byte = (unsigned char)(0.0);
	    else if (help > 255.0)
		byte = (unsigned char)(255.0);
	    else
		byte = (unsigned char)(help);
	    fwrite (&byte, sizeof(unsigned char), 1, file);
	}

    /* close file */
    fclose(file);

    printf("... SUCCESS");
}


/* ------------------------------------------------------------------------- */

void read_barron_data
(
                            /************************************************/
    char  *filename,        /* in : file name                               */
    float **u,              /* in : x-component of vector data              */
    float **v,              /* in : y-component of vector data              */
    int   nx,               /* in : size in x-direction                     */
    int   ny,               /* in : size in y-direction                     */
    int   bx,               /* in : boundary in x-direction                 */
    int   by                /* in : boundary in y-direction                 */
                            /************************************************/
)

/* reads barron file */

{
                           /*************************************************/
    FILE *file;            /* file pointer                                  */
    float help;            /* tmp variable                                  */
    float **uref;          /* tmp array                                     */
    float **vref;          /* tmp array                                     */
    int i,j;               /* loop variabeles                               */
    int nxref_and_offsetx; /* size in x-direction with crop offset          */
    int nyref_and_offsety; /* size in y-direction with crop offset          */
    int nxref_no_offsetx;  /* size in x-direction without crop offset       */
    int nyref_no_offsety;  /* size in y-direction without crop offset       */
    int offsetx;           /* crop offset in x-direction                    */
    int offsety;           /* crop offset in y-direction                    */
                           /*************************************************/


    printf("\n Trying to read barron %s ...",filename);

    /* try to open file */
    file = fopen(filename,"r");

    /* if file not found */
    if (file==NULL)
    {
	printf("... FAILED");
	printf("\n\n PROGRAM ABORTED !!! \n\n");
	exit(0);
    }

    /* read header */
    fread (&help, sizeof(float), 1, file);
    nxref_and_offsetx  = (int) help;
    fread (&help, sizeof(float), 1, file);
    nyref_and_offsety  = (int) help;
    fread (&help, sizeof(float), 1, file);
    nxref_no_offsetx  = (int) help;
    fread (&help, sizeof(float), 1, file);
    nyref_no_offsety  = (int) help;
    fread (&help, sizeof(float), 1, file);
    offsetx = (int) help;
    fread (&help, sizeof(float), 1, file);
    offsety = (int) help;

    /* compare dimensions for consistency */
    if ((nx!=nxref_no_offsetx)||(ny!=nyref_no_offsety))
    {
	printf("... WRONG DIMENSIONS");
	printf("\n\n PROGRAM ABORTED !!! \n\n");
	exit(0);
    }

    /* allocate memory for tmp array */
    ALLOC_MATRIX(2,nxref_and_offsetx,nyref_and_offsety,&uref,&vref);

    /* read barron data */
    for (j = 0; j < nyref_and_offsety; j++)
	for (i = 0; i < nxref_and_offsetx; i++)
	{
	    fread(&help, sizeof(float), 1, file);
	    uref[i][j] = (float) help;
	    fread(&help, sizeof(float), 1, file);
	    vref[i][j] = (float) help;
	}


    /* copy data without cropped border */
    for (i=bx; i<nx+bx; i++)
	for (j=by; j<ny+by; j++)
	{
	    u[i][j] = (float) uref[i-bx+offsetx][j-by+offsety];
	    v[i][j] = (float) vref[i-bx+offsetx][j-by+offsety];
	}

    /* free memory for tmp array */
    FREE_MATRIX(2,nxref_and_offsetx,nyref_and_offsety,uref,vref);


    /* close file */
    fclose(file);

    printf("... SUCCESS");
}

/* ------------------------------------------------------------------------- */

void write_barron_data
(
                            /************************************************/
    char  *filename,        /* in : file name                               */
    float **u,              /* in : x-component of vector data              */
    float **v,              /* in : y-component of vector data              */
    int   nx,               /* in : size in x-direction                     */
    int   ny,               /* in : size in y-direction                     */
    int   bx,               /* in : boundary in x-direction                 */
    int   by                /* in : boundary in y-direction                 */
                            /************************************************/
)

/* writes barron file */

{
                  /**********************************************************/
    FILE *file;   /* file pointer                                           */
    float help;   /* tmp variable                                           */
    int i,j;      /* loop variables                                         */
    int offset;   /* border size to crop (set fixed to 0)                   */
                  /**********************************************************/

    printf("\n Trying to write barron file %s ...",filename);

    /* try to open file */
    file = fopen(filename,"w");

    /* if file not found */
    if (file==NULL)
    {
	printf("... FAILED");
	printf("\n\n PROGRAM ABORTED !!! \n\n");
	exit(0);
    }

    /* write header */
    help = nx;
    fwrite (&help, 4, 1, file);
    help = ny;
    fwrite (&help, 4, 1, file);
    offset=0;
    help = nx - 2 * offset;
    fwrite (&help, 4, 1, file);
    help = ny - 2 * offset;
    fwrite (&help, 4, 1, file);
    help = offset;
    fwrite (&help, 4, 1, file);
    fwrite (&help, 4, 1, file);

    /* write data */
    for (j=by; j<ny+by; j++)
	for (i=bx; i<nx+bx; i++)
	{
	    help = (float)u[i][j];
	    fwrite (&help, 4, 1, file);

	    help = (float)v[i][j];
	    fwrite (&help, 4, 1, file);
	}

    /* close file */
    fclose(file);

    printf("... SUCCESS");
}


/* ------------------------------------------------------------------------- */


void store_displacements_to_image
(
                            /************************************************/
    float  **u,             /* in  : x component of flow field              */
    float  **v,             /* in  : y component of flow field              */
    float  **f1,            /* in  : 1st image                              */
    int    nx,              /* in  : size in x-direction                    */
    int    ny,              /* in  : size in y-direction                    */
    int    bx,              /* in  : boundary size in x-direction           */
    int    by,              /* in  : boundary size in y-direction           */
    float  min_length,      /* in  : treshold for maximum flow              */
    float  max_length,      /* in  : treshold for minimum flow              */
    float  ***p6            /* out : RGB image                              */
                            /************************************************/
)

/* creates image from displacement field*/

{
                  /**********************************************************/
    int   i,j;    /* loop variable                                          */
    int   R,G,B;  /* RGB values                                             */
    float factor; /* scaling factor                                         */
    float tmp;    /* aux variable                                           */
                  /**********************************************************/

    /* compute scaling factor */
    factor=1.0/max_length;

    /* original grey value image with color flow field */
    for(i=bx; i < nx+bx;i++)
      for(j=by; j < ny+by;j++)
    	{
    	  /* left side : original image */
    	  p6[i][j][0]=byte_range((int)u[i][j]);
    	  p6[i][j][1]=byte_range((int)v[i][j]);
    	  p6[i][j][2]=byte_range((int)f1[i][j]);
    	}
}





/* ------------------------------------------------------------------------- */


float* convert_displacements_to_image(
                            /************************************************/
    float  **u,             /* in  : x component of flow field              */
    float  **v,             /* in  : y component of flow field              */
    int    nx,              /* in  : size in x-direction                    */
    int    ny,              /* in  : size in y-direction                    */
    int    bx,              /* in  : boundary size in x-direction           */
    int    by,              /* in  : boundary size in y-direction           */
    float  min_length,      /* in  : treshold for maximum flow              */
    float  max_length      /* in  : treshold for minimum flow              */
                            /************************************************/
)

/* creates image from displacement field*/

{
                  /**********************************************************/
    int   i,j;    /* loop variable                                          */
    int   R,G,B;  /* RGB values                                             */
    float factor; /* scaling factor                                         */
    float tmp;    /* aux variable                                           */
                  /**********************************************************/

    /* compute scaling factor */
    factor=1.0/max_length;
    float * img= new float[3*nx*ny];

    /* original grey value image with color flow field */
    for(i=bx; i < nx+bx;i++)
      for(j=by; j < ny+by;j++)
    	{
    	  /* compute RGB value for current flow vector */
    	  vector_to_RGB(u[i][j]*factor,v[i][j]*factor,&R,&G,&B);
        img[(j-by)*nx*3+(i-bx)*3] = byte_range(R)/255.0;
        img[(j-by)*nx*3+(i-bx)*3+1] = byte_range(G)/255.0;
        img[(j-by)*nx*3+(i-bx)*3+2] = byte_range(B)/255.0;
    	}
    return img;
}


/* ------------------------------------------------------------------------- */


void copy_refocus_image
(
                            /************************************************/
    float  **f1,            /* in  : 1st image                              */
    float  ***refocus,      /* in  : refocus image                          */
    int    nx,              /* in  : size in x-direction                    */
    int    ny,              /* in  : size in y-direction                    */
    int    bx,              /* in  : boundary size in x-direction           */
    int    by,              /* in  : boundary size in y-direction           */
    float  min_length,      /* in  : treshold for maximum flow              */
    float  max_length,      /* in  : treshold for minimum flow              */
    float  ***p6            /* out : RGB image                              */
                            /************************************************/
)

/* creates image from displacement field*/

{
                  /**********************************************************/
    int   i,j;    /* loop variable                                          */
    int   R,G,B;  /* RGB values                                             */
    float factor; /* scaling factor                                         */
    float tmp;    /* aux variable                                           */
                  /**********************************************************/

    /* compute scaling factor */
    factor=1.0/max_length;

    /* original grey value image with color flow field */
    tmp = 0;
    for(i=bx; i < nx+bx;i++)
      for(j=by; j < ny+by;j++)
	{
	  /* left side : original image */
	  p6[i][j][0]=byte_range((int)f1[i][j]);
	  p6[i][j][1]=byte_range((int)f1[i][j]);
	  p6[i][j][2]=byte_range((int)f1[i][j]);

	  /* compute RGB value for current flow vector */
	  //vector_to_RGB(u[i][j]*factor,v[i][j]*factor,&R,&G,&B);
    if(tmp<refocus[i][j][0])
      tmp = refocus[i][j][0];
    R= (int)(refocus[i][j][0]*255);
    G= (int)(refocus[i][j][1]*255);
    B= (int)(refocus[i][j][2]*255);


	  /* right side : flow field */
	  p6[i+nx][j][0]=(int)byte_range(R);
	  p6[i+nx][j][1]=(int)byte_range(G);
	  p6[i+nx][j][2]=(int)byte_range(B);
	}
  //printf("max channel 0 %f\n",tmp);

}


/* ------------------------------------------------------------------------- */


void convert_displacements_to_image
(
                            /************************************************/
    float  **u,             /* in  : x component of flow field              */
    float  **v,             /* in  : y component of flow field              */
    float  **f1,            /* in  : 1st image                              */
    int    nx,              /* in  : size in x-direction                    */
    int    ny,              /* in  : size in y-direction                    */
    int    bx,              /* in  : boundary size in x-direction           */
    int    by,              /* in  : boundary size in y-direction           */
    float  min_length,      /* in  : treshold for maximum flow              */
    float  max_length,      /* in  : treshold for minimum flow              */
    float  ***p6            /* out : RGB image                              */
                            /************************************************/
)

/* creates image from displacement field*/

{
                  /**********************************************************/
    int   i,j;    /* loop variable                                          */
    int   R,G,B;  /* RGB values                                             */
    float factor; /* scaling factor                                         */
    float tmp;    /* aux variable                                           */
                  /**********************************************************/

    /* compute scaling factor */
    factor=1.0/max_length;

    /* original grey value image with color flow field */
    for(i=bx; i < nx+bx;i++)
      for(j=by; j < ny+by;j++)
	{
	  /* left side : original image */
	  p6[i][j][0]=byte_range((int)f1[i][j]);
	  p6[i][j][1]=byte_range((int)f1[i][j]);
	  p6[i][j][2]=byte_range((int)f1[i][j]);

	  /* compute RGB value for current flow vector */
	  vector_to_RGB(u[i][j]*factor,v[i][j]*factor,&R,&G,&B);

	  /* right side : flow field */
	  p6[i+nx][j][0]=(int)byte_range(R);
	  p6[i+nx][j][1]=(int)byte_range(G);
	  p6[i+nx][j][2]=(int)byte_range(B);
	}
}


/* ------------------------------------------------------------------------- */

inline void print_console_line
(
)
{
printf("\n -----------------------------------------------------------------");
}
/* ------------------------------------------------------------------------- */
inline void print_console_header
(
char *header_text
)
{
//printf("\n");
print_console_line();
printf("\n %s",header_text);
print_console_line();
}

/* ------------------------------------------------------------------------- */
#endif
