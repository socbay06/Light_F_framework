#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <fstream>
#include <sys/time.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <iostream>
#include "debug.h"
#include "./hdf5_lib/hdf5_tools.h"
#include "./png_lib/libPNG.h"


/*---------------------------------------------------------------------------*/
/* include own libraries */
#include "./misc_lib/alloc_mem_linear.c"
#include "./misc_lib/alloc_mem_linear_mult.c"
#include "./misc_lib/io_lib.c"
#include "./misc_lib/bounds_lib.c"
#include "./misc_lib/matrix_lib.c"
#include "./misc_lib/conv_lib.c"
#include "./misc_lib/funct_lib.c"
#include "./misc_lib/of_lib.c"
#include "./misc_lib/string_lib.c"
#include "./misc_lib/mg_trans_lib.c"

#include "./horn_schunck/horn_schunck_warp_view.c"

/*---------------------------------------------------------------------------*/
#define MODE_MANUAL 0
#define MODE_INTERACTIVE 1

//define arguments variable.
char in_h5_file[160]; // input hdf5 file path.
char dataset[160]; // dataset name in h5 file.
char out_h5_file[160];
char config_file[160];


/*---------------------------------------------------------------------------*/
struct variation_flow_params {
   float g_pixel_ratio;
	float g_max_disp;
	float h_sigma;
	float h_alpha;
	float h_beta;
	int g_type;
	float h_omega;
	int h_warp_levels;
	int h_iter;
	float h_warp_eta;
} hs_param; //horn_schnuck paramater

struct light_field_params {
	//data related params.
	float * lf;
   size_t W;
	size_t H;
	size_t S;
	size_t T;
	size_t C;
} lf;

struct application_params {
	int mode = MODE_MANUAL;
	int png_out = 1;
	int view_all = 0;
	int view_t = -1;
	int view_s = -1;
} app_param;


/****************************************************/

char g_image1[80];       /* name of 1st image                                */
char g_image2[80];       /* name of 2nd image                                */
char g_ref[80];          /* name of problem ground truth file                */
char lf_h5[80];          /* name of lightfield h5 file                       */
char lf_dataset[160];    /* name of dataset                                  */
                         /*                                                  */
float  **g_f1;           /* 1st image                                        */
float  **g_f2;           /* 2nd image                                        */
float  **g_f1_s;         /* smoothed 1st image                               */
float  **g_f2_s;         /* smoothed 2nd image                               */
float  **g_f2_bw;        /* motion compensated 2nd image                     */
float  **g_u;            /* displacement component in x-direction            */
float  **g_v;            /* displacement component in y-direction            */
float  **g_u_ref;        /* problem ground truth in x-direction              */
float  **g_v_ref;        /* problem ground truth in y-direction              */
                         /*                                                  */
float  ***g_p6;          /* colour array to draw image                       */
float  ***g_disp;          /* colour array to draw image                       */
                           /*                                                  */
int g_nx,g_ny;           /* size of both frames in x- and y-drection         */
int g_bx,g_by;           /* size of image boundary                           */
float g_hx,g_hy;         /* pixel size of grid                               */
float * uvbuffer;       /* Buffer for writing out u/v to h5 file */

/*                                                  */
float g_ns;
float g_nt;
float g_m_beta;

float g_e_sigma;         /* standard deviation of the Gaussian presmoothing  */
                         /*                                                  */
int g_m_type;                /* type of constancy assumption             */
float g_m_alpha;         /* weight of the smoothness term                    */
int g_n_iter;            /* number of iterations                             */
float g_n_omega;         /* SOR overrelaxtion parameter                      */
float g_n_warp_eta;      /* warping rescaling between levels                 */
int g_n_warp_levels;     /* number of warping levels                         */
                         /*                                                  */
float g_g_pixel_ratio;   /* pixel zoom factor                                */
float g_g_max_disp;      /* maxium displacement magnitude                    */
                         /*                                                  */
                         /*                                                  */
int g_active_param;      /* active parameter                                 */
int g_direct_compute;    /* flag for direct computation                      */
                         /*                                                  */
                         /*                                                  */
float g_aae;             /* average angualer error                           */
float g_al2e;            /* average l2 norm error                            */
float g_ref_density;     /* density of problem ground truth                  */
float g_est_density;     /* density of estimation                            */
                         /*                                                  */
                         /*                                                  */
long g_position;         /* position in stream for reading input             */
                         /*                                                  */
                         /*                                                  */
GLubyte *g_pixels;       /* pixel aray for Open-GL                           */
                         /****************************************************/



/*---------------------------------------------------------------------------*/


void drawGlutScene_from_image
(
        /****************************************************/
        float   ***p6,   /* in  : RGB image                                  */
        GLubyte *pixels, /* use : display array                              */
        float magnify,   /* in  : scaling factor                             */
        int nx,          /* in  : size in x-direction                        */
        int ny,          /* in  : size in y-direction                        */
        int bx,          /* in  : boundary size in x-direction               */
        int by           /* in  : boundary size in y-direction               */
                         /****************************************************/
)

/* visualises image with Open-GL */

{

        int odd;  /* flag for odd image size */
        int counter; /* pixel index counter */
        int i,j;  /* loop variable */


        /* check if image size is odd */
        if (nx%4==0) odd=0;
        else odd=1;

        /* set pixel counter zero */
        counter=0;

        /* prepare Open-GL */
        glViewport(0, 0, nx, ny);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0, nx, 0, ny, -1.0, 1.0);
        glMatrixMode(GL_MODELVIEW);
        glDisable(GL_DITHER);
        glPixelZoom((GLfloat)magnify,(GLfloat)magnify);

        /* draw pixels in pixel array */
        for(i=by; i < ny+by; i++)
        {
                for(j=bx; j < nx+bx; j++)
                {
                        pixels[counter++]=(GLubyte)(p6[j][ny+2*by-1-i][0]);
                        pixels[counter++]=(GLubyte)(p6[j][ny+2*by-1-i][1]);
                        pixels[counter++]=(GLubyte)(p6[j][ny+2*by-1-i][2]);
                }
                if (odd==1) counter+=2;
        }

        /* draw pixels in draw buffer */
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearColor(0.0, 0.0, 0.0, 1.0);
        glRasterPos3f(0, 0, 0.0);
        glDrawPixels(nx,ny,GL_RGB,GL_UNSIGNED_BYTE,pixels);

        /* swap draw and display buffer */
        glutSwapBuffers();

        return;
}

/*------------------------------------------------*/

float* lf_extract_view(light_field_params data,  int view_t, int view_s){
	size_t W,H,S,T,C;
	W=data.W;
	H=data.H;
	T=data.T;
	S=data.S;
	C=data.C;
	float* view = new float[W*H*C];
   int center = floor((S-1)/2.0);
   for(int j = 0; j<H; j++) {
   	for(int i =0; i<W; i++) {
     		for(int c=0; c<C; c++) {
         	view[j*W*C+i*C+c]=data.lf[c*W*H*S*T + i*H*S*T + j*S*T + (center+view_s)*T + (center+view_t)];
         }
      }
   }
	return view;
}

float* lf_extract_view_gray(light_field_params lf, int view_t, int view_s){
	float* view = lf_extract_view(lf,view_t,view_s);
	float* view_gray = convertRGBtoGray(lf.W,lf.H,view);
	free(view);
	return view_gray;
}


/*---------------------------------------------------------------------------*/


void print_constancy_type(){
        if(g_m_type==TYPE_GRAY)
                printf("\n \t constancy: GRAY ");
        else if(g_m_type==TYPE_HESSIAN)
                printf("\n \t constancy: HESSIAN ");
        else
                printf("\n \t constancy: GRADIENT ");
}


void showParams()
{
        print_console_header("Model");

        printf("\n Horn/Schunck Approach");

        print_console_header("Model Parameters");

        if (g_active_param==1)
                printf("\n (a) (Alpha)                     %4.6lf",(double)g_m_alpha);
        else
                printf("\n (a)  Alpha                      %4.6lf",(double)g_m_alpha);


        print_console_header("Preprocessing");

        if (g_active_param==0)
                printf("\n (p) (Sigma)                     %4.3lf",(double)g_e_sigma);
        else
                printf("\n (p)  Sigma                      %4.3lf",(double)g_e_sigma);

        print_console_header("Numerical Parameters");

        if (g_active_param==2)
                printf("\n (i) (Iterations)                %d",g_n_iter);
        else
                printf("\n (i)  Iterations                 %d",g_n_iter);

        if (g_active_param==5)
                printf("\n (o) (Omega SOR)                 %4.3lf",(double)g_n_omega);
        else
                printf("\n (o)  Omega SOR                  %4.3lf",(double)g_n_omega);

        if (g_active_param==6)
                printf("\n (w) (Warping Levels)            %d",g_n_warp_levels);
        else
                printf("\n (w)  Warping Levels             %d",g_n_warp_levels);

        if (g_active_param==7)
                printf("\n (e) (Warping Eta)               %4.3lf",(double)g_n_warp_eta);
        else
                printf("\n (e)  Warping Eta                %4.3lf",(double)g_n_warp_eta);
        if (g_active_param==8)
                printf("\n (t) (Constancy type)               %d",g_m_type);
        else
                printf("\n (t)  Constancy type                %d",g_m_type);
        print_constancy_type();
        print_console_header("Visualisation Parameters");

        if (g_active_param==9)
                printf("\n (b) (beta view smoothness)               %4.6lf",g_m_beta);
        else
                printf("\n (b)  beta view smoothness               %4.6lf",g_m_beta);

        if (g_active_param==4)
                printf("\n (g) (Maximum Shown Disp)        %4.3lf",(double)g_g_max_disp);
        else
                printf("\n (g)  Maximum Shown Disp         %4.3lf",(double)g_g_max_disp);

        if (g_direct_compute==1)
                printf("\n (,)  Direct Computation         ON");
        else
                printf("\n (,)  Direct Computation         OFF");

        print_console_header("Quality Measures");

        if(g_u_ref!=NULL)
                printf("\n      Average Angular Error      %4.10lf",(double)g_aae);

        print_console_line();
        fflush(stdout);
        return;
}


/*---------------------------------------------------------------------------*/
//
//
// void handleComputeDisplacements()
// {
//
//         /* presmooth both frames */
//         presmooth_2d(g_f1,g_f1_s,g_nx,g_ny,g_bx,g_by,g_hx,g_hy,g_e_sigma,g_e_sigma);
//         presmooth_2d(g_f2,g_f2_s,g_nx,g_ny,g_bx,g_by,g_hx,g_hy,g_e_sigma,g_e_sigma);
//
//         printf("after presmooth_2d\n");
//
//         /* initialise displacement vector field with zero */
//         set_matrix_2d(g_u,g_nx+2*g_bx,g_ny+2*g_by,0,0,(float)0.0);
//         set_matrix_2d(g_v,g_nx+2*g_bx,g_ny+2*g_by,0,0,(float)0.0);
//
//         printf("after init matrix\n");
//
//         /* set boundaries zero */
//         set_bounds_2d(g_u,g_nx,g_ny,g_bx,g_by,(float)0.0);
//         set_bounds_2d(g_v,g_nx,g_ny,g_bx,g_by,(float)0.0);
//         set_bounds_2d(g_f1_s,g_nx,g_ny,g_bx,g_by,(float)0.0);
//         set_bounds_2d(g_f2_s,g_nx,g_ny,g_bx,g_by,(float)0.0);
//
//         printf("after set bounds\n");
//         /* compute displacement vector field */
//         // HORN_SCHUNCK_MAIN(g_f1_s,g_f2_s,g_u,g_v,g_nx,g_ny,g_bx,g_by,
//         //                   g_hx, g_hy, g_m_type, g_m_alpha, g_n_iter, g_n_omega,
//         //                   g_n_warp_eta, g_n_warp_levels);
//         /* compute displacement vector field */
//         HORN_SCHUNCK_VIEW(g_f1_s,g_f2_s,g_u,g_v,g_nx,g_ny,g_bx,g_by,
//                           g_hx, g_hy, g_ns, g_nt, g_m_beta, g_m_alpha, g_n_iter);
//         printf("after schunk\n");
//
//         /* if ground truth available compute average angular error */
//         if(g_u_ref!=NULL)
//                 calculate_errors_2d(g_u_ref,g_v_ref,g_u,g_v,g_nx,g_ny,g_bx,g_by,
//                                     &g_aae,&g_al2e,&g_ref_density,&g_est_density);
//
//         return;
// }


/*---------------------------------------------------------------------------*/


void handleComputeDisplacements()
{

        /* presmooth both frames */
        presmooth_2d(g_f1,g_f1_s,g_nx,g_ny,g_bx,g_by,g_hx,g_hy,g_e_sigma,g_e_sigma);
        presmooth_2d(g_f2,g_f2_s,g_nx,g_ny,g_bx,g_by,g_hx,g_hy,g_e_sigma,g_e_sigma);

        printf("after presmooth_2d\n");

        /* initialise displacement vector field with zero */
        set_matrix_2d(g_u,g_nx+2*g_bx,g_ny+2*g_by,0,0,(float)0.0);
        set_matrix_2d(g_v,g_nx+2*g_bx,g_ny+2*g_by,0,0,(float)0.0);

        printf("after init matrix\n");

        /* set boundaries zero */
        set_bounds_2d(g_u,g_nx,g_ny,g_bx,g_by,(float)0.0);
        set_bounds_2d(g_v,g_nx,g_ny,g_bx,g_by,(float)0.0);
        set_bounds_2d(g_f1_s,g_nx,g_ny,g_bx,g_by,(float)0.0);
        set_bounds_2d(g_f2_s,g_nx,g_ny,g_bx,g_by,(float)0.0);

        printf("after set bounds\n");
        /* compute displacement vector field */
        HORN_SCHUNCK_MAIN(g_f1_s,g_f2_s,g_u,g_v,g_nx,g_ny,g_bx,g_by,
                          g_hx, g_hy, g_m_type, g_ns,g_nt,g_m_beta, g_m_alpha, g_n_iter, g_n_omega,
                          g_n_warp_eta, g_n_warp_levels);
        //HORN_SCHUNCK_VIEW(g_f1_s,g_f2_s,g_u,g_v,g_nx,g_ny,g_bx,g_by,
        //                   g_hx, g_hy,g_ns,g_nt,g_m_beta, g_m_alpha, g_n_iter);  //from horn_
        printf("after schunk\n");

        /* if ground truth available compute average angular error */
        if(g_u_ref!=NULL)
                calculate_errors_2d(g_u_ref,g_v_ref,g_u,g_v,g_nx,g_ny,g_bx,g_by,
                                    &g_aae,&g_al2e,&g_ref_density,&g_est_density);

        return;
}


/*---------------------------------------------------------------------------*/

void handleDraw()
{

        /* create image */
        convert_displacements_to_image(g_u,g_v,g_f1_s,g_nx,g_ny,g_bx,g_by,
                                       (float)0.0,g_g_max_disp,g_p6);

        /* store displacement */
        store_displacements_to_image(g_u,g_v,g_f1_s,g_nx,g_ny,g_bx,g_by,
                                     (float)0.0,g_g_max_disp,g_disp);

        /* draw displacement field field in pixel array */
        drawGlutScene_from_image(g_p6,g_pixels,g_g_pixel_ratio,
                                 2*g_nx,g_ny,g_bx,g_by);

        /* draw pixel array on screen */
        glutPostRedisplay();
}


/*---------------------------------------------------------------------------*/

void handleKeyboardspecial(int key, int x, int y)
{
        /* keyboard handler */

        switch(key) {
        case GLUT_KEY_F1:
                /* call help menu */
                print_console_line();
                printf("\n\n F1     ........this help\n");
                printf(" F2     ........write out displacement field in colour code \n");
                printf(" F8     ........write out motion compensated 2nd frame \n");
                printf(" ESC    ........program termination\n");
                printf(" p      ........select presmoothing parameter\n");
                printf(" a      ........select smoothness weight\n");
                printf(" i      ........select number of iterations\n");
                printf(" o      ........select SOR overrelaxation parameter\n");
                printf(" w      ........select number of warping levels\n");
                printf(" e      ........select warping rescaling factor\n");
                printf(" g      ........select max displacement magnitude\n");
                printf(" up     ........increase active parameter\n");
                printf(" down   ........decrease active parameter\n");
                printf(" ,      ........direct computation on/off \n");
                printf(" .      ........compute displacement field \n");
                print_console_line();
                fflush(stdout);
                break;

        case GLUT_KEY_F2:
                /* write out displacement field */
                if (g_u_ref!=NULL)
                {
                        int i,j;
                        for (i=g_bx; i<g_nx+g_bx; i++)
                                for (j=g_by; j<g_ny+g_by; j++)
                                {
                                        if((g_u_ref[i][j]==100.0)&&(g_v_ref[i][j]==100.0))
                                        {
                                                g_u[i][j]=100.0;
                                                g_v[i][j]=100.0;
                                        }
                                }
                }
                write_ppm_blank_header("out.ppm",g_nx,g_ny);
                write_ppm_data("out.ppm",g_p6,g_nx,g_ny,g_bx+g_nx,g_by);
                write_ppm_blank_header("displacement.ppm",g_nx,g_ny);
                write_ppm_data("displacement.ppm",g_disp,g_nx,g_ny,g_bx,g_by);
                break;

        case GLUT_KEY_F8:
                /* write out motion compensated 2nd frame */
                backward_registration(g_f1,g_f2,g_f2_bw,g_u,g_v,
                                      g_nx,g_ny,g_bx,g_by,g_hx,g_hy);
                write_pgm_blank_header("frame2_bw.pgm",g_nx,g_ny);
                write_pgm_data("frame2_bw.pgm",g_f2_bw,g_nx,g_ny,g_bx,g_by);
                break;

        case GLUT_KEY_DOWN:
                /* decrease sigma */
                if (g_active_param==0)
                {
                        g_e_sigma=g_e_sigma/1.1;
                        if (g_e_sigma<0.3)
                                g_e_sigma=0.0;
                        if (g_direct_compute==1) handleComputeDisplacements(); break;
                }
                /* decrease alpha */
                if (g_active_param==1)
                {
                        g_m_alpha/=sqrt(sqrt(sqrt(sqrt(10))));
                        if (g_m_alpha<0.01) g_m_alpha=0.01;
                        if (g_direct_compute==1) handleComputeDisplacements(); break;
                }
                /* decrease beta */
                if (g_active_param==9)
                {
                        g_m_beta/=sqrt(sqrt(sqrt(sqrt(10))));
                        if (g_m_beta<0.01) g_m_beta=0.01;
                        if (g_direct_compute==1) handleComputeDisplacements(); break;
                }
                /* decrease number of iterations */
                if (g_active_param==2)
                {
                        g_n_iter-=100;
                        if (g_n_iter<0)
                                g_n_iter=0;
                        if (g_direct_compute==1) handleComputeDisplacements(); break;
                }
                /* decrease displacement field magnitude (visualsiation) */
                if (g_active_param==4)
                {
                        if (g_g_max_disp>=0.01)
                                g_g_max_disp=g_g_max_disp/1.1;
                        break;
                }
                /* decrease omega */
                if (g_active_param==5)
                {
                        g_n_omega=g_n_omega-0.01;
                        if (g_n_omega<0.01)
                                g_n_omega=0.01;
                        if (g_direct_compute==1) handleComputeDisplacements(); break;
                }
                /* decrease number of warping levels */
                if (g_active_param==6)
                {
                        g_n_warp_levels-=1;
                        if (g_n_warp_levels<0)
                                g_n_warp_levels=0;
                        if (g_direct_compute==1) handleComputeDisplacements(); break;
                }
                /* decrease type */
                if (g_active_param==8)
                {
                        g_m_type-=1;
                        if (g_m_type<0)
                                g_m_type=0;
                        if (g_direct_compute==1) handleComputeDisplacements(); break;
                }
                /* decrease eta */
                if (g_active_param==7)
                {
                        g_n_warp_eta-=0.05;
                        if (g_n_warp_eta<0.05)
                                g_n_warp_eta=0.05;
                        if (g_direct_compute==1) handleComputeDisplacements(); break;
                }

        case GLUT_KEY_UP:
                /* increase sigma */
                if (g_active_param==0)
                {
                        if (g_e_sigma<0.30)
                                g_e_sigma=0.3;
                        else
                                g_e_sigma=g_e_sigma*1.1;
                        if (g_direct_compute==1) handleComputeDisplacements(); break;
                }
                /* increase alpha */
                if (g_active_param==1)
                {
                        g_m_alpha*=sqrt(sqrt(sqrt(sqrt(10))));
                        if (g_direct_compute==1) handleComputeDisplacements(); break;
                }
                /* increase beta */
                if (g_active_param==9)
                {
                        g_m_beta*=sqrt(sqrt(sqrt(sqrt(10))));
                        if (g_direct_compute==1) handleComputeDisplacements(); break;
                }
                /* increase number of iterations */
                if (g_active_param==2)
                {
                        g_n_iter+=100;
                        if (g_direct_compute==1) handleComputeDisplacements(); break;
                }
                /* increase displacement field magnitude (visualisation) */
                if (g_active_param==4)
                {
                        g_g_max_disp=g_g_max_disp*1.1;
                        break;
                }
                /* increase omega */
                if (g_active_param==5)
                {
                        g_n_omega+=0.01;
                        if (g_n_omega>1.99)
                                g_n_omega=1.99;
                        if (g_direct_compute==1) handleComputeDisplacements(); break;
                }
                /* decrease number of warping levels */
                if (g_active_param==6)
                {
                        g_n_warp_levels+=1;
                        if (g_direct_compute==1) handleComputeDisplacements(); break;
                }
                /* increase eta */
                if (g_active_param==7)
                {
                        g_n_warp_eta+=0.05;
                        if (g_n_warp_eta>0.95)
                                g_n_warp_eta=0.95;
                        if (g_direct_compute==1) handleComputeDisplacements(); break;
                }
                /* increase type */
                if (g_active_param==8)
                {
                        g_m_type+=1;
                        if (g_m_type>10)
                                g_m_type=10;
                        if (g_direct_compute==1) handleComputeDisplacements(); break;
                }
        default:
                printf("\nUnknown key pressed (Code %d).",key);

        }

        /* show new parameters and displacement field */
        showParams();
        handleDraw();
        return;
}


/*---------------------------------------------------------------------------*/

void handleKeyboard(unsigned char key, int x, int y)
{
        /* keyboard handler */
        switch(key) {
        case 'p':
                g_active_param=0; break;
        case 't':
                g_active_param=8; break;
        case 'a':
                g_active_param=1; break;
        case 'i':
                g_active_param=2; break;
        case 'g':
                g_active_param=4; break;
        case 'o':
                g_active_param=5; break;
        case 'w':
                g_active_param=6; break;
        case 'e':
                g_active_param=7; break;
        case 'b':
                g_active_param=9; break;
        case 44: //,
                if (g_direct_compute==1) {g_direct_compute=0; break; }
                if (g_direct_compute==0) {g_direct_compute=1; break; }
        case 46: //.
                handleComputeDisplacements(); break;
        case 27:
                exit(1);
        default:
                printf("\nUnknown key pressed (Code %d).",key);
        }

        /* show new parameters and displacement field */
        showParams();
        handleDraw();
        return;
}

/*---------------------------------------------------------------------------*/

void handleComputeNothing()
{
}

/*---------------------------------------------------------------------------*/

void handleMouse(int button, int state, int cx, int cy)
{
        printf("\nNo mouse handler yet.");
}



/** ESTIMATE FLOW FOR TWO IMAGES AND WRITE TO HDF5 FILE **/
// img1: first image;
// img2: second image;
// filename: image file to write out estimated flow.
// h5outfile: hdf5 file to write out estimated flow.
// h5dset: dataset inside hdf5 file for storing estimated flow.
void calculate_flow_view_write_img(float* img1, float* img2, char* filename, char* h5outfile,char* h5dset){
	 //copy view to image pair.
    int dims[10];

    for(int j=0; j<g_ny; j++) {
    	for(int i=0; i<g_nx; i++) {
      	g_f1[i+g_bx][j+g_by]= img1[j*g_nx+i]*255;
         g_f2[i+g_bx][j+g_by]= img2[j*g_nx+i]*255;
      }
    }
	fprintf(stdout,"*************************************************\n");
   fprintf(stdout, " Start to calculate displacement for %s \n", h5dset);
	fprintf(stdout,"*************************************************\n");

   handleComputeDisplacements();

	//printf("Finish calculate flowfield, write to image \n");
	float* flow = convert_displacements_to_image(g_u,g_v,g_nx,g_ny,g_bx,g_by,(float)0.0,g_g_max_disp);
   //printf("Finish convert to rgb %d %d %d %d \n", g_nx, g_ny, g_bx, g_by);
	writeImageRGB(filename,g_nx,g_ny,flow,filename);
 	free(flow);

   //write to hdf5 output
   dims[0] = 2; //u or v
   dims[1] = g_nx; // width;
   dims[2] = g_ny; //height;
   //copy flow to buffer*

   for(int j=0; j<g_ny; j++) {
   	for(int i=0; i<g_nx; i++) {
      	uvbuffer[ i*g_ny + j ] = g_u[g_bx+i][g_by+j];
         uvbuffer[ g_ny*g_nx + i*g_ny +j ]=g_v[g_bx+i][g_by+j];
      }
	}

   hdf5_write_data_simple(h5outfile,h5dset,dims, 3, uvbuffer);
}


/*-------------- FUll flow field estimation ------*/
// Generate displacement for all posible views
void full_flow_view_estimation(light_field_params data, char* out_file){
	float normfactor;
   // RGB image from three views;
   char name1[256];
   char name2[256];
	int S,T,W,H,C;
   // grayscale image from three view;
   float * img_c;//view center in grayscale
   float * img_x;//secondary view in grayscale

	S = data.S;
	T = data.T;
	W = data.W;
	H = data.H;
	C = data.C;

   // the center view index.
   int center = floor((S-1)/2.0);
   int vpatch = center -1; //Ignore view at the border

   //prepare view center
   img_c = lf_extract_view_gray(data,0,0);
   writeImageGray("view_c.png",W,H,img_c,"Image center");
	fprintf(stdout, " EXTRACT VIEW: view center to view_c.png\n");
   // let do with other view
   for(int vs=-vpatch; vs<=vpatch; vs++) {
   	for(int vt=-vpatch; vt<=vpatch; vt++) {
      	if(vs==0 && vt==0)
         	continue;
			//Calculate the directional terms.
         normfactor = 1.0/sqrt(vs*vs+vt*vt);
         g_ns = normfactor*vs;
         g_nt = normfactor*vt;

			img_x = lf_extract_view_gray(data,vt,vs);
         sprintf(name1,"view_%d_%d.png",vt,vs);
         sprintf(name2,"Image at view (T,S) (%d,%d)",vt,vs);
         writeImageGray(name1,W,H,img_x,name2);
			fprintf(stdout," EXTRACT VIEW: view (T,S)  %d,%d to %s \n",vt,vs,name1);
         sprintf(name2,"/flow/v_%d_%d",vt,vs);
         sprintf(name1,"flow_%d_%d.png",vt,vs);
         calculate_flow_view_write_img(img_c,img_x,name1,out_file,name2);
		}
	}
   free(img_x);
   free(img_c);
}



void handle_function_manual(){
	fprintf(stdout,"Handle Function in Manual mode\n");
	full_flow_view_estimation(lf,out_h5_file);
}

void handle_function_interactive( int* argc, char** argv){
	float* view_c;
	float* view_o;
	int s;
	int t;
	int center;


	fprintf(stdout,"Handle Functions in Interactive mode\n");

	view_c = lf_extract_view_gray(lf,0,0);
	writeImageGray("view_c.png",lf.W,lf.H,view_c,"Image center");
	fprintf(stdout,"Extracted view center to view_c.png\n");

	s = app_param.view_s;
	t = app_param.view_t;
	center = floor((lf.S-1)/2.0);

	if(abs(s)>center)
		s = floor(center/2.0);
	if(abs(t)>center)
		t = floor(center/2.0);

	view_o = lf_extract_view_gray(lf,t,s);//other view.
   writeImageGray("view_o.png",lf.W,lf.H,view_o,"Image secondary view");
	fprintf(stdout,"Extracted view %d,%d to view_o.png\n",t,s);

   for(int j=0; j<g_ny; j++) {
   	for(int i=0; i<g_nx; i++) {
      	g_f1[i+g_bx][j+g_by]= view_c[j*g_nx+i]*255;
         g_f2[i+g_bx][j+g_by]= view_o[j*g_nx+i]*255;
      }
   }
	float	normfactor = 1.0/sqrt(s*s+t*t);
   g_ns = normfactor*s;
   g_nt = normfactor*t;

// open OpenGL window */
   glutInit(argc, argv);
   glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

   glutInitWindowSize((int)round(2*g_nx*g_g_pixel_ratio),
                           (int)round(g_ny*g_g_pixel_ratio));
   glutCreateWindow("CORRESPONDENCE PROBLEMS - VISUALISATION FRONTEND");

// register handle routines
   glutDisplayFunc(handleDraw);
   glutIdleFunc(handleComputeNothing);
   glutKeyboardFunc(handleKeyboard);
   glutSpecialFunc(handleKeyboardspecial);
   glutMouseFunc(handleMouse);
//
// main
   handleComputeDisplacements();
   handleDraw();
   showParams();
   glutMainLoop();
}




void lf_load_configuration(char* config_file){

	//At first, we should load the default configuration.
	hs_param.g_pixel_ratio = 1;
	hs_param.g_max_disp = 2.0;
	hs_param.h_sigma = 1.2;
	hs_param.h_alpha = 100;
	hs_param.h_beta = 100;
	hs_param.g_type = 0;
	hs_param.h_omega = 1.95;
	hs_param.h_warp_levels = 6;
	hs_param.h_iter = 200;
	hs_param.h_warp_eta = 0.5;


	// start to read from config file.
	if(strlen(config_file)!=0){
		fprintf(stdout,"Reading configuration file %s \n",config_file);
 		std::ifstream file(config_file);
   	std::string str;
		std::string params;
		std::size_t pos;
		std::string  values;
   	while (std::getline(file, str))
   	{
			fprintf(stdout,"line: %s \n",str.c_str());
		   trim(str);
			if(str.empty())
				continue;
			params = str.substr(0,1);
			if(params.compare("#")==0)
				continue;
			pos = str.find("#");
			if(pos != std::string::npos)
				str = str.substr(0,pos);
			pos = str.find(";");
			if(pos != std::string::npos)
				str = str.substr(0,pos);
			if(str.empty())
				continue;
			pos = str.find("=");
			if(pos == std::string::npos)
				continue;
			params = str.substr(0,pos);
			trim(params);
			values = str.substr(pos+1,str.length()-(pos+1));
			trim(values);
			fprintf(stdout,"Got param %s values %s \n",params.c_str(), values.c_str());
	//read parameter for variation flow algorithm.
			if(params.compare("g_pixel_ratio")==0)
				hs_param.g_pixel_ratio = atoi(values.c_str());
			if(params.compare("g_max_disp")==0)
				hs_param.g_max_disp = atof(values.c_str());
			if(params.compare("h_sigma")==0)
				hs_param.h_sigma = atof(values.c_str());
			if(params.compare("h_alpha")==0)
				hs_param.h_alpha = atof(values.c_str());
			if(params.compare("h_beta")==0)
				hs_param.h_beta = atof(values.c_str());
			if(params.compare("g_type")==0)
				hs_param.g_type = atoi(values.c_str());
			if(params.compare("h_omega")==0)
				hs_param.h_omega = atof(values.c_str());
 			if(params.compare("h_warp_levels")==0)
				hs_param.h_warp_levels = atof(values.c_str());
			if(params.compare("h_warp_eta")==0)
				hs_param.h_warp_eta = atof(values.c_str());
			if(params.compare("h_iter")==0)
				hs_param.h_iter = atoi(values.c_str());

	// Applicaton level params.
			if(params.compare("view_t")==0)
				app_param.view_t = atoi(values.c_str());
 			if(params.compare("view_s")==0)
				app_param.view_s = atoi(values.c_str());
			if(params.compare("view_all")==0)
				app_param.view_all = atoi(values.c_str());
			if(params.compare("p_mode")==0)
				app_param.mode = atoi(values.c_str());
 			if(params.compare("png_out")==0)
				app_param.png_out = atoi(values.c_str());
   	}
	}

	//transfer configuration to original variables.
   g_g_pixel_ratio = hs_param.g_pixel_ratio;
   g_g_max_disp=hs_param.g_max_disp;
   g_direct_compute=0;
   g_e_sigma=hs_param.h_sigma;
   g_m_alpha=hs_param.h_alpha;
   g_m_beta = hs_param.h_beta;
   g_m_type=hs_param.g_type;
   g_active_param=1000;
   g_u_ref=NULL;
   g_v_ref=NULL;
   g_n_omega=hs_param.h_omega;
   g_n_warp_levels=hs_param.h_warp_levels;
   g_n_warp_eta=hs_param.h_warp_eta;
   g_n_iter = hs_param.h_iter; // number of iteration

}

void print_usage(){
        fprintf(stdout,"USAGE: ./lf -m [0,1] -i <file> -d <string> [-o <file>] [-h] [-c <file>]\n");
        fprintf(stdout,"\t m: mode. 0 - manual | 1 - interactive\n");
        fprintf(stdout,"\t i: input hdf5 file\n");
        fprintf(stdout,"\t d: dataset name\n");
        fprintf(stdout,"\t o: output hdf5 file\n");
        fprintf(stdout,"\t h: this help\n");
        fprintf(stdout,"\t c: config file\n");
}


int main (int argc, char* argv[]){
   char option;
	int mode_arg; // mode from command line arg;
	//set initial params.
	strcpy(in_h5_file,"");
	strcpy(dataset,"");
	mode_arg=-1;
	strcpy(out_h5_file,"");
	strcpy(config_file,"");
        while ((option = getopt(argc, argv,"i:c:ho:m:d:")) != -1) {
                switch (option) {
                case 'h':
                        print_usage();
                        exit(0);
                        break;
                case 'i': //Input HDF5 File.
                        strcpy(in_h5_file,optarg);
                        break;
                case 'd': //dataset name
                        strcpy(dataset,optarg);
                        break;
                case 'm':
                        mode_arg = atoi(optarg);
                        break;
                case 'o':
                        strcpy(out_h5_file,optarg);
                        break;
                case 'c':
                        strcpy(config_file,optarg);
                        break;
                default:
                        print_usage();
                        exit(0);
                }
        }

	//verify parameters.
	//input hdf5file.
	if(strlen(in_h5_file)==0){
		fprintf(stdout," Please provide input h5 file\n");
		print_usage();
		exit(0);
	}

	//check dataset name.
	if(strlen(dataset) == 0){
		fprintf(stdout," Please provide dataset name\n");
		print_usage();
		exit(0);
	}
	fprintf(stdout," Working with H5 file %s and dataset %s \n",in_h5_file, dataset);

	//check ouptut file
	if(strlen(out_h5_file)==0){
		strcpy(out_h5_file,"output.h5");
	}
	fprintf(stdout," Output file is %s \n",out_h5_file);


	//read configuration file
	lf_load_configuration(config_file);

	if(mode_arg!=-1)
		app_param.mode=mode_arg;
	if(app_param.mode == MODE_MANUAL)
		fprintf(stdout," Working mode is set to  manual");
	else{
		app_param.mode = MODE_INTERACTIVE;
		fprintf(stdout," Working mode is set to interactive");
	}
	//Check input HDF5 file.
   fprintf(stdout," Loading HDF5 file: %s\n",in_h5_file);
   fprintf(stdout," Reading LF from dataset %s\n",dataset);
   lf.lf = hdf5_read_lf(in_h5_file,dataset,lf.W,lf.H,lf.S,lf.T,lf.C);
   fprintf(stdout," .... LF dimension: %dx%dx%dx%dx%d\n",lf.W,lf.H,lf.S,lf.T,lf.C);

/* ---------- set boundary and grid size ----------------------------------- */

        g_bx=1;
        g_by=1;

        g_hx=1;
        g_hy=1;


	//setup image size for algorithms.
   g_nx=lf.W;
   g_ny=lf.H;

   /* ---------- memory allocation ------------------------------------ */

   ALLOC_MATRIX(2, g_nx+2*g_bx, g_ny+2*g_by, &g_f1,   &g_f2);
   ALLOC_MATRIX(2, g_nx+2*g_bx, g_ny+2*g_by, &g_f1_s, &g_f2_s);
   ALLOC_MATRIX(1, g_nx+2*g_bx, g_ny+2*g_by, &g_f2_bw);
   ALLOC_MATRIX(2, g_nx+2*g_bx, g_ny+2*g_by, &g_u,    &g_v);
   ALLOC_CUBIX(1, 2*g_nx+2*g_bx, g_ny+2*g_by, 3, &g_p6);
   ALLOC_CUBIX(1, g_nx+2*g_bx, g_ny+2*g_by, 3, &g_disp);

   uvbuffer = new float[g_nx*g_ny*2];
   g_pixels = new GLubyte[2*(g_nx+1)*g_ny*3*sizeof(GLubyte)];

	/*------ Extract view ------ */


	// the center view index.


	if(app_param.mode == MODE_MANUAL)
		handle_function_manual();
	else
		handle_function_interactive(&argc,argv);



/* ---------- free memory -------------------------------------------------- */

        FREE_MATRIX (2, g_nx+2*g_bx, g_ny+2*g_by, g_f1,   g_f2);
        FREE_MATRIX (2, g_nx+2*g_bx, g_ny+2*g_by, g_f1_s, g_f2_s);
        FREE_MATRIX (1, g_nx+2*g_bx, g_ny+2*g_by, g_f2_bw);
        FREE_MATRIX (2, g_nx+2*g_bx, g_ny+2*g_by, g_u,    g_v);
        //
        // printf("free matrix  \n");
        //
        // //TODO handle ground truth.
        // // if (argc==5)
        // //         FREE_MATRIX(2,   g_nx+2*g_bx, g_ny+2*g_by, g_u_ref, g_v_ref);
        //
        FREE_CUBIX (1, 2*g_nx+2*g_bx, g_ny+2*g_by, 3, g_p6);
        FREE_CUBIX (1, g_nx+2*g_bx, g_ny+2*g_by, 3, g_disp);


        free(g_pixels);
        free(uvbuffer);

        printf("\n\n\n");


        return(0);

}



int main0 (int argc, char* argv[])
{
        int dims[10];
        char modeParams[3];
        char out_file[80];
        int mode =0;
        size_t W,H,S,T,C;                        /*                                                  */
        float* data;      /* lightfield data read from hdf5 file */

/* ---------- set boundary and grid size ----------------------------------- */

        g_bx=1;
        g_by=1;

        g_hx=1;
        g_hy=sqrt(3)/2;

/* ---------- read in arguments -------------------------------------------- */


        if (argc <4) {
                fprintf(stderr, "Not sufficient arguments: \n");
                fprintf(stderr, "<APP> <mode> <LF H5 in file> <LF in dataset> <LF H5 out file>\n");
                fprintf(stderr, " mode: 0-interactive 1-commandline");
                exit(1);
        }
        mode = MODE_INTERACTIVE;
        strcpy(modeParams, argv[1]);
        mode = atoi(modeParams);
        strcpy(lf_h5,argv[2]);
        strcpy(lf_dataset,argv[3]);
        strcpy(out_file,argv[4]);
        /*TODO Group truth*/
        //strcpy(g_ref,argv[4]);


/* set default parameters */
        g_g_pixel_ratio = 1;
        g_g_max_disp=2.0;
        g_direct_compute=0;
        g_e_sigma=1.2;
        g_m_alpha=100;
        g_m_beta = 100;
        g_m_type=0;
        g_active_param=1000;
        g_u_ref=NULL;
        g_v_ref=NULL;
        g_n_omega=1.95;
        g_n_warp_levels=6;
        g_n_warp_eta=0.5;
        g_n_iter = 200; // number of iteration


/* ---------- read in information of light field ------------------- */

        data = hdf5_read_lf(lf_h5, lf_dataset,W,H,S,T,C);

        //exit(1);
        g_nx=W;
        g_ny=H;

        printf(" after read lf \n");

        /* ---------- memory allocation ------------------------------------ */

        ALLOC_MATRIX(2, g_nx+2*g_bx, g_ny+2*g_by, &g_f1,   &g_f2);
        ALLOC_MATRIX(2, g_nx+2*g_bx, g_ny+2*g_by, &g_f1_s, &g_f2_s);
        ALLOC_MATRIX(1, g_nx+2*g_bx, g_ny+2*g_by, &g_f2_bw);
        ALLOC_MATRIX(2, g_nx+2*g_bx, g_ny+2*g_by, &g_u,    &g_v);
//          printf(" after alloc 1 \n");
        //if (argc==5)
        //        ALLOC_MATRIX(2, g_nx+2*g_bx, g_ny+2*g_by, &g_u_ref,   &g_v_ref);

        ALLOC_CUBIX(1, 2*g_nx+2*g_bx, g_ny+2*g_by, 3, &g_p6);
//        printf(" after alloc 2 \n");
        ALLOC_CUBIX(1, g_nx+2*g_bx, g_ny+2*g_by, 3, &g_disp);

//        printf(" after alloc 3 \n");
        //g_pixels = (GLubyte *) malloc (2*(g_nx+1)*g_ny*3*sizeof(GLubyte));
        uvbuffer = new float[g_nx*g_ny*2];
        g_pixels = new GLubyte[2*(g_nx+1)*g_ny*3*sizeof(GLubyte)];

        /* full view flow estimation */
        //full_flow_view_estimation(data, W,H,S,T,C, out_file);

        /* extract two view */
        // RGB image from three views;
        float * view00 = new float[W*H*3];
        float * view03 = new float[W*H*3];
        float * view0_3 = new float[W*H*3];
        float * view30 = new float[W*H*3];
        float * view_30 = new float[W*H*3];
        // grayscale image from three view;
        float * img_c;
        float * img_l;
        float * img_r;
        float * img_u;
        float * img_d;
        // the center view index.
        int center = floor((S-1)/2.0);
        for(int j = 0; j<H; j++) {
                for(int i =0; i<W; i++) {
                        for(int c=0; c<3; c++) {
                                view00[j*W*3+i*3+c]=data[c*W*H*S*T + i*H*S*T + j*S*T + (center)*T + center];
                                view03[j*W*3+i*3+c]=data[c*W*H*S*T + i*H*S*T + j*S*T + (center+2)*T + center];
                                view0_3[j*W*3+i*3+c]=data[c*W*H*S*T + i*H*S*T + j*S*T + (center+2)*T + center];
                                view30[j*W*3+i*3+c]=data[c*W*H*S*T + i*H*S*T + j*S*T + (center+3)*T + center];
                                view_30[j*W*3+i*3+c]=data[c*W*H*S*T + i*H*S*T + j*S*T + (center+4)*T + center];
                        }
                }
        }
        img_c = convertRGBtoGray(W,H,view00);
        img_l = convertRGBtoGray(W,H,view_30);
        img_r = convertRGBtoGray(W,H,view30);
        img_d = convertRGBtoGray(W,H,view0_3);
        img_u = convertRGBtoGray(W,H,view03);
        free(view00);
        free(view03);
        free(view0_3);
        free(view30);
        free(view_30);
        /* Write gray image  */
        writeImageGray("view_c.png",W,H,img_c,"Image center");
        writeImageGray("view_l.png",W,H,img_l,"Image Left");
        writeImageGray("view_r.png",W,H,img_r,"Image Right");
        writeImageGray("view_u.png",W,H,img_u,"Image Up");
        writeImageGray("view_d.png",W,H,img_d,"Image Down");
//
//
//
//
// /* ---------- read in image pair ------------------------------------------- */
//         //
//         calculate_flow_write_img(img_c,img_l,"flow_left.png",out_file,"/flow/left");
//         calculate_flow_write_img(img_c,img_r,"flow_right.png",out_file,"/flow/right");
//         calculate_flow_write_img(img_c,img_u,"flow_up.png",out_file,"/flow/up");
//         calculate_flow_write_img(img_c,img_d,"flow_down.png",out_file,"/flow/down");
/* ---------- read in ground truth displacement field ---------------------- */

        //if (argc==5)
        //        read_barron_data(g_ref,g_u_ref,g_v_ref,g_nx,g_ny,g_bx,g_by);

//
/* ---------- M A I N   L O O P -------------------------------------------- */

        for(int j=0; j<g_ny; j++) {
                for(int i=0; i<g_nx; i++) {
                        g_f1[i+g_bx][j+g_by]= img_c[j*g_nx+i]*255;
                        g_f2[i+g_bx][j+g_by]= img_u[j*g_nx+i]*255;
                }
        }
        g_ns=1.0;
        g_nt=0; //1.0/sqrt(2);

// open OpenGL window */
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

        glutInitWindowSize((int)round(2*g_nx*g_g_pixel_ratio),
                           (int)round(g_ny*g_g_pixel_ratio));
        glutCreateWindow("CORRESPONDENCE PROBLEMS - VISUALISATION FRONTEND");

// register handle routines
        glutDisplayFunc(handleDraw);
        glutIdleFunc(handleComputeNothing);
        glutKeyboardFunc(handleKeyboard);
        glutSpecialFunc(handleKeyboardspecial);
        glutMouseFunc(handleMouse);
//
// main
        handleComputeDisplacements();
        handleDraw();
        showParams();
        glutMainLoop();




/* ---------- free memory -------------------------------------------------- */

        FREE_MATRIX (2, g_nx+2*g_bx, g_ny+2*g_by, g_f1,   g_f2);
        FREE_MATRIX (2, g_nx+2*g_bx, g_ny+2*g_by, g_f1_s, g_f2_s);
        FREE_MATRIX (1, g_nx+2*g_bx, g_ny+2*g_by, g_f2_bw);
        FREE_MATRIX (2, g_nx+2*g_bx, g_ny+2*g_by, g_u,    g_v);
        //
        // printf("free matrix  \n");
        //
        // //TODO handle ground truth.
        // // if (argc==5)
        // //         FREE_MATRIX(2,   g_nx+2*g_bx, g_ny+2*g_by, g_u_ref, g_v_ref);
        //
        FREE_CUBIX (1, 2*g_nx+2*g_bx, g_ny+2*g_by, 3, g_p6);
        FREE_CUBIX (1, g_nx+2*g_bx, g_ny+2*g_by, 3, g_disp);


        free(g_pixels);
        free(uvbuffer);

        printf("\n\n\n");


        return(0);
}
