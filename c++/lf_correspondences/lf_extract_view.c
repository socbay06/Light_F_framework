
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <fstream>
#include <sys/time.h>

#include <iostream>
#include "debug.h"
#include "./hdf5_lib/hdf5_tools.h"
#include "./png_lib/libPNG.h"


/*---------------------------------------------------------------------------*/
#define MODE_MANUAL 0
#define MODE_INTERACTIVE 1

//define arguments variable.
char in_h5_file[160]; // input hdf5 file path.
char dataset[160]; // dataset name in h5 file.
char out_h5_file[160];
char config_file[160];


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
        int view_s = 0;
        int view_t = 0;
        bool full_view = false;
        bool color = true;
} app_param;


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


void print_usage(){
        fprintf(stdout,"USAGE: ./lf -i <file> -d <string> [-o <folder>] [-h] [-c,g] [-s int] [-t int] [-a]\n");
        fprintf(stdout,"\t i: input hdf5 file\n");
        fprintf(stdout,"\t d: dataset name\n");
        fprintf(stdout,"\t o: output hdf5 file\n");
        fprintf(stdout,"\t h: this help\n");
        fprintf(stdout,"\t c: Output color image \n");
        fprintf(stdout,"\t g: Output gray image file\n");
        fprintf(stdout,"\t a: Extract all views \n");
        fprintf(stdout,"\t s: view s\n");
        fprintf(stdout,"\t t: view t\n");
}


int main (int argc, char* argv[]){
        char option;
        int mode_arg; // mode from command line arg;
        int view_s;
        int view_t;
        bool color; // extract color image;
        bool extract_all;
        float * img;
        char filename[128];
        //set initial params.
        strcpy(in_h5_file,"");
        strcpy(dataset,"");
        strcpy(out_h5_file,"");
        strcpy(config_file,"");
        view_s = 0;
        view_t = 0;
        extract_all = false;
        color =true;

        while ((option = getopt(argc, argv,"i:d:o:gchas:t:")) != -1) {
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
                case 's':
                        view_s = atoi(optarg);
                        break;
                case 't':
                        view_t = atoi(optarg);
                        break;
                case 'o':
                        strcpy(out_h5_file,optarg);
                        break;
                case 'c':
                        color = true;
                        break;
                case 'g':
                        color = false;
                        break;
                case 'a':
                        extract_all = true;
                        break;
                default:
                        print_usage();
                        exit(0);
                }
        }

        //verify parameters.
        //input hdf5file.
        if(strlen(in_h5_file)==0) {
                fprintf(stdout," Please provide input h5 file\n");
                print_usage();
                exit(0);
        }

        //check dataset name.
        if(strlen(dataset) == 0) {
                fprintf(stdout," Please provide dataset name\n");
                print_usage();
                exit(0);
        }
        fprintf(stdout," Working with H5 file %s and dataset %s \n",in_h5_file, dataset);

        //TODO check ouptut folder
        // ....
        // ....

        //Check input HDF5 file.
        fprintf(stdout," Loading HDF5 file: %s\n",in_h5_file);
        fprintf(stdout," Reading LF from dataset %s\n",dataset);
        lf.lf = hdf5_read_lf(in_h5_file,dataset,lf.W,lf.H,lf.S,lf.T,lf.C);
        fprintf(stdout," .... LF dimension: %dx%dx%dx%dx%d\n",lf.W,lf.H,lf.S,lf.T,lf.C);

        if(extract_all) {
                fprintf(stdout,"Extracting all views TxS = %d X %d ... \n", lf.T, lf.S);
                for(int j =0; j<lf.T; j++) {
                        for(int i=0; i<lf.S; i++) {
                                view_t = j;
                                view_s = i;
                                fprintf(stdout,"\t\t View (t,s) = (%d,%d) .. \n", view_t,view_s);
                                sprintf(filename,"view_%d_%d.png",view_t,view_s);
                                if(color) {
                                        img = lf_extract_view(lf,view_t,view_s);
                                        writeImageRGB(filename,lf.W,lf.H,img,filename);
                                }else{
                                        img = lf_extract_view_gray(lf,view_t,view_s);
                                        writeImageGray(filename,lf.W,lf.H,img,filename);
                                }
                                free(img);
                        }
                }
        }else{
                fprintf(stdout,"Extracting view (t,s) = (%d,%d) .... \n",view_t,view_s);
                sprintf(filename,"view_%d_%d.png",view_t,view_s);
                if(color) {
                        img = lf_extract_view(lf,view_t,view_s);
                        writeImageRGB(filename,lf.W,lf.H,img,filename);
                }else{
                        img = lf_extract_view_gray(lf,view_t,view_s);
                        writeImageGray(filename,lf.W,lf.H,img,filename);
                }
                free(img);
        }

        return(0);

}
