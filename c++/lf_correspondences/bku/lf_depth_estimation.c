#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <iostream>
#include "debug.h"
#include "hdf5_tools.h"
#include "libPNG.h"


/*---------------------------------------------------------------------------*/
/* include own libraries */

#include "alloc_mem_linear_mult.c"
#include "io_lib.c"
#include "bounds_lib.c"
#include "matrix_lib.c"
#include "conv_lib.c"
#include "funct_lib.c"
#include "of_lib.c"
#include "horn_schunck_warp_view.c"

/*---------------------------------------------------------------------------*/

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




void calculate_flow_view_write_img(float* img1, float* img2, char* filename, char* h5outfile,char* h5dset){

        //copy view to image pair.
        int dims[10];

        for(int j=0; j<g_ny; j++) {
                for(int i=0; i<g_nx; i++) {
                        g_f1[i+g_bx][j+g_by]= img1[j*g_nx+i]*255;
                        g_f2[i+g_bx][j+g_by]= img2[j*g_nx+i]*255;
                }
        }

        printf(" Start to calculate displacement for %s \n", h5dset);
        handleComputeDisplacements();


        //printf("Finish calculate flowfield, write to image \n");
        float* flow = convert_displacements_to_image(g_u,g_v,g_nx,g_ny,g_bx,g_by,(float)0.0,g_g_max_disp);
        //printf("Finish convert to rgb %d %d %d %d \n", g_nx, g_ny, g_bx, g_by);
        writeImageRGB(filename,g_nx,g_ny,flow,filename);
        free(flow);

        //write to output
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
        //printf("write ok  \n");
}


void calculate_flow_write_img(float* img1, float* img2, char* filename, char* h5outfile,char* h5dset){

        //copy view to image pair.
        int dims[10];

        for(int j=0; j<g_ny; j++) {
                for(int i=0; i<g_nx; i++) {
                        g_f1[i+g_bx][j+g_by]= img1[j*g_nx+i]*255;
                        g_f2[i+g_bx][j+g_by]= img2[j*g_nx+i]*255;
                }
        }

        printf(" Start to calculate displacement for %s \n", filename);
        handleComputeDisplacements();


        //printf("Finish calculate flowfield, write to image \n");
        float* flow = convert_displacements_to_image(g_u,g_v,g_nx,g_ny,g_bx,g_by,(float)0.0,g_g_max_disp);
        //printf("Finish convert to rgb %d %d %d %d \n", g_nx, g_ny, g_bx, g_by);
        writeImageRGB(filename,g_nx,g_ny,flow,filename);
        free(flow);

        //write to output
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
        //printf("write ok  \n");
}



/*-------------- FUll flow field estimation ------*/
void full_flow_view_estimation(float * data, int W, int H, int S, int T, int C, char* out_file){
        float normfactor;
        // RGB image from three views;
        float * viewC = new float[W*H*3];
        float * viewX = new float[W*H*3];
        char name1[256];
        char name2[256];

        // grayscale image from three view;
        float * img_c;
        float * img_x;

        // the center view index.
        int center = floor((S-1)/2.0);
        int vpatch = center -2;


        //prepare view center
        for(int j = 0; j<H; j++) {
                for(int i =0; i<W; i++) {
                        for(int c=0; c<3; c++) {
                                viewC[j*W*3+i*3+c]=data[c*W*H*S*T + i*H*S*T + j*S*T + (center)*T + center];
                        }
                }
        }
        img_c = convertRGBtoGray(W,H,viewC);
        writeImageGray("view_c.png",W,H,img_c,"Image center");

        // let do with other view
        for(int vs=-vpatch; vs<vpatch; vs++)
                for(int vt=-vpatch; vt<vpatch; vt++) {
                        if(vs==0 && vt==0)
                          continue;
                        //extracts view from light field.
                        for(int j = 0; j<H; j++) {
                                for(int i =0; i<W; i++) {
                                        for(int c=0; c<3; c++) {
                                                viewX[j*W*3+i*3+c]=data[c*W*H*S*T + i*H*S*T + j*S*T + (center+vs)*T + center+vt];
                                        }
                                }
                        }
                        normfactor = 1.0/sqrt(vs*vs+vt*vt);
                        g_ns = normfactor*vs;
                        g_nt = normfactor*vt;
                        img_x = convertRGBtoGray(W,H,viewX);
                        sprintf(name1,"view_%d_%d.png",vs,vt);
                        sprintf(name2,"Image at view (S,T) (%d,%d)",vs,vt);
                        writeImageGray(name1,W,H,img_x,name2);
                        sprintf(name2,"/flow/v_%d_%d",vs,vt);
                        sprintf(name1,"flow_%d_%d.png",vs,vt);
                        calculate_flow_view_write_img(img_c,img_x,name1,out_file,name2);
                }
        free(viewC);
        free(viewX);
        free(img_x);
        free(img_c);
}



/*---------------------------------------------------------------------------*/

int main (int argc, char* argv[])
{
        int dims[10];
        char out_file[80];
        size_t W,H,S,T,C;                        /*                                                  */
        float* data;      /* lightfield data read from hdf5 file */

/* ---------- set boundary and grid size ----------------------------------- */

        g_bx=1;
        g_by=1;

        g_hx=1;
        g_hy=sqrt(3)/2;

/* ---------- read in arguments -------------------------------------------- */


        if (argc <3) {
                fprintf(stderr, "Not sufficient arguments: \n");
                fprintf(stderr, "<APP> <LF H5 in file> <LF in dataset> <LF H5 out file>\n");
                exit(1);
        }
        strcpy(lf_h5,argv[1]);
        strcpy(lf_dataset,argv[2]);
        strcpy(out_file,argv[3]);
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

        g_n_iter = 800; // number of iteration




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
        g_nt=0;//1.0/sqrt(2);

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
