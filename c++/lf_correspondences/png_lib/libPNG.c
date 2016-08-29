// LibPNG example
// A.Greensted
// http://www.labbookpages.co.uk

// Version 2.0
// With some minor corrections to Mandlebrot code (thanks to Jan-Oliver)

// Version 1.0 - Initial release
#include "libPNG.h"

// int main(int argc, char *argv[])
// {
// 	// Make sure that the output filename argument has been provided
// 	if (argc != 2) {
// 		fprintf(stderr, "Please specify output file\n");
// 		return 1;
// 	}
//
// 	// Specify an output image size
// 	int width = 500;
// 	int height = 300;
//
// 	// Create a test image - in this case a Mandelbrot Set fractal
// 	// The output is a 1D array of floats, length: width * height
// 	printf("Creating Image\n");
// 	float *buffer = createMandelbrotImage(width, height, -0.802, -0.177, 0.011, 110);
// 	if (buffer == NULL) {
// 		return 1;
// 	}
//
// 	// Save the image to a PNG file
// 	// The 'title' string is stored as part of the PNG file
// 	printf("Saving PNG\n");
// 	int result = writeImage(argv[1], width, height, buffer, "This is my test image");
//
// 	// Free up the memorty used to store the image
// 	free(buffer);
//
// 	return result;
// }

inline void setRGB(png_byte *ptr, float val)
{
	int v = (int)(val * 767);
	if (v < 0) v = 0;
	if (v > 767) v = 767;
	int offset = v % 256;

	if (v<256) {
		ptr[0] = 0; ptr[1] = 0; ptr[2] = offset;
	}
	else if (v<512) {
		ptr[0] = 0; ptr[1] = offset; ptr[2] = 255-offset;
	}
	else {
		ptr[0] = offset; ptr[1] = 255-offset; ptr[2] = 0;
	}
}




int writeImageGray(char* filename, int width, int height, float *buffer, char* title)
{
	int code = 0;
	FILE *fp = NULL;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_bytep row = NULL;

	// Open file for writing (binary mode)
	fp = fopen(filename, "wb");
	if (fp == NULL) {
		fprintf(stderr, "Could not open file %s for writing\n", filename);
		code = 1;
		goto finalise;
	}

	// Initialize write structure
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (png_ptr == NULL) {
		fprintf(stderr, "Could not allocate write struct\n");
		code = 1;
		goto finalise;
	}

	// Initialize info structure
	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL) {
		fprintf(stderr, "Could not allocate info struct\n");
		code = 1;
		goto finalise;
	}

	// Setup Exception handling
	if (setjmp(png_jmpbuf(png_ptr))) {
		fprintf(stderr, "Error during png creation\n");
		code = 1;
		goto finalise;
	}

	png_init_io(png_ptr, fp);

	// Write header (8 bit colour depth)
	png_set_IHDR(png_ptr, info_ptr, width, height,
			8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
			PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	// Set title
	if (title != NULL) {
		png_text title_text;
		title_text.compression = PNG_TEXT_COMPRESSION_NONE;
		title_text.key = "Title";
		title_text.text = title;
		png_set_text(png_ptr, info_ptr, &title_text, 1);
	}

	png_write_info(png_ptr, info_ptr);

	// Allocate memory for one row (3 bytes per pixel - RGB)
	//row = (png_bytep) malloc(width * sizeof(png_byte));
	row = new png_byte[width * sizeof(png_byte)];
	// Write image data
	int x, y;
	for (y=0 ; y<height ; y++) {
		for (x=0 ; x<width ; x++) {
			row[x]=(int)(buffer[y*width+x]*255);
		}
		png_write_row(png_ptr, row);
	}

	// End write
	png_write_end(png_ptr, NULL);

	finalise:
	if (fp != NULL) fclose(fp);
	if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	if (row != NULL) free(row);

	return code;
}




int writeImageRGB(char* filename, int width, int height, float *buffer, char* title)
{
	int code = 0;
	FILE *fp = NULL;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_bytep row = NULL;
	float maxval = 65535;

	// Open file for writing (binary mode)
	fp = fopen(filename, "wb");
	if (fp == NULL) {
		fprintf(stderr, "Could not open file %s for writing\n", filename);
		code = 1;
		goto finalise;
	}

	// Initialize write structure
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (png_ptr == NULL) {
		fprintf(stderr, "Could not allocate write struct\n");
		code = 1;
		goto finalise;
	}

	// Initialize info structure
	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL) {
		fprintf(stderr, "Could not allocate info struct\n");
		code = 1;
		goto finalise;
	}

	// Setup Exception handling
	if (setjmp(png_jmpbuf(png_ptr))) {
		fprintf(stderr, "Error during png creation\n");
		code = 1;
		goto finalise;
	}

	png_init_io(png_ptr, fp);

	// Write header (8 bit colour depth)
	png_set_IHDR(png_ptr, info_ptr, width, height,
			8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
			PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	// Set title
	if (title != NULL) {
		png_text title_text;
		title_text.compression = PNG_TEXT_COMPRESSION_NONE;
		title_text.key = "Title";
		title_text.text = title;
		png_set_text(png_ptr, info_ptr, &title_text, 1);
	}

	png_write_info(png_ptr, info_ptr);

	// Allocate memory for one row (3 bytes per pixel - RGB)
	//row = (png_bytep) malloc(3 * width * sizeof(png_byte));
	row = new png_byte[3 * width * sizeof(png_byte)];
	// Write image data
	int x, y;
	float to;
	to=0;
	for (y=0 ; y<height ; y++) {
		for (x=0 ; x<width ; x++) {
			//setRGB(&(row[x*3]), buffer[y*width + x]);
			if(buffer[y*width+x+2]>to)
				to=buffer[y*width+x+2];
			row[x*3] = (png_byte)1.0*buffer[y*width*3+3*x+0]*255.0;
			row[x*3+1] = (png_byte)1.0*buffer[y*width*3+3*x+1]*255.0;
			row[x*3+2] = (png_byte)1.0* buffer[y*width*3+3*x+2]*255.0;
		}
		png_write_row(png_ptr, row);
	}

	printf("max channel 1 %f \n",to);

	// End write
	png_write_end(png_ptr, NULL);

	finalise:
	if (fp != NULL) fclose(fp);
	if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	if (row != NULL) free(row);

	return code;
}



float *createMandelbrotImage(int width, int height, float xS, float yS, float rad, int maxIteration)
{
	//float *buffer = (float *) malloc(width * height * sizeof(float));
	float *buffer = new float[width * height * sizeof(float)];
	if (buffer == NULL) {
		fprintf(stderr, "Could not create image buffer\n");
		return NULL;
	}

	// Create Mandelbrot set image

	int xPos, yPos;
	float minMu = maxIteration;
	float maxMu = 0;

	for (yPos=0 ; yPos<height ; yPos++)
	{
		float yP = (yS-rad) + (2.0f*rad/height)*yPos;

		for (xPos=0 ; xPos<width ; xPos++)
		{
			float xP = (xS-rad) + (2.0f*rad/width)*xPos;

			int iteration = 0;
			float x = 0;
			float y = 0;

			while (x*x + y*y <= 4 && iteration < maxIteration)
			{
				float tmp = x*x - y*y + xP;
				y = 2*x*y + yP;
				x = tmp;
				iteration++;
			}

			if (iteration < maxIteration) {
				float modZ = sqrt(x*x + y*y);
				float mu = iteration - (log(log(modZ))) / log(2);
				if (mu > maxMu) maxMu = mu;
				if (mu < minMu) minMu = mu;
				buffer[yPos * width + xPos] = mu;
			}
			else {
				buffer[yPos * width + xPos] = 0;
			}
		}
	}

	// Scale buffer values between 0 and 1
	int count = width * height;
	while (count) {
		count --;
		buffer[count] = (buffer[count] - minMu) / (maxMu - minMu);
	}

	return buffer;
}


float * convertRGBtoGray(int width, int height,float* origin){
	float * gray_img = new float[width * height];
	float R,G,B,est;
	for(int j =0; j<height; j++)
		for(int i=0; i<width; i++){
			R = sRGB_to_linear(origin[j*width*3+i*3]);
			G = sRGB_to_linear(origin[j*width*3+i*3+1]);
			B = sRGB_to_linear(origin[j*width*3+i*3+2]);
			est = 0.2126 * R + 0.7152 * G + 0.0722*B;
			gray_img[j*width+i] = linear_to_sRGB(est);
		}

	return gray_img;
}



float sRGB_to_linear(float x){
	if(x<0.04045) return x/12.92;
	return pow((x+0.055)/1.055,2.4);
}
/* Apply gramma correction*/
float linear_to_sRGB(float y){
	if(y<=0.0031308) return 19.92*y;
	return 1.055*pow(y,1/2.4)-0.055;
}
