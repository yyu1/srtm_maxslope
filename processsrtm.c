//This program will process the global 1arcsec SRTM mosaic, which is ~1.2T in size

//The program utilizes OpenMP to process the data on multiple threads

//The SRTM 1arcsec mosaic is 16-bit integer type with dimension of 1296000x475200
//Upper left corner coordinate is -180 Longitude, 76 Latitude  (lat/lon projection)

//LDAN nodes 4-8 has 256GB of memory, if we divide into 10 blocks, each block is
//47520 rows of data, which is divisible by 3, so we can do our 3x3 averaging.
//1296000x47520 gives  123,171,840,000 : ~ 123GB, so we can allocate 125GB heap block
//To read in the chunk of input data, and allocate a 13GB (123/9) block for writing output.

//For each 123GB block, which is 47520 rows.  LDAN has Sandy bridge processor with 2x8 core = 16 cores
//We can further divide the 47520 rows into 20 sub-blocks, and process them in parallel using OpenMP

#include <stdio.h>
#include <string.h>
#include <time.h>
#include "maxdegreeslope.h"

//Image size and blocking settings
#define IMAGEXDIM 1296000ULL
#define IMAGEYDIM 475200ULL
#define BLOCKS 10ULL
#define SUBBLOCKS 20ULL



int main() {
	
	char infile[100] = "/nobackupp6/nexprojects/CMS-ALOS/srtm1sec_mosaic/global_srtm3_aster_mosaic_1sec.int";
	char outfile[100] = "/nobackupp6/nexprojects/CMS-ALOS/srtm1sec_mosaic/global_srtm3_aster_maxdegreeslope_3sec.byt";
	char logfile[100] = "/nobackupp6/nexprojects/CMS-ALOS/srtm1sec_mosaic/maxslope_3sec.log";

	//Open log file
	FILE* log_file_ptr;
	log_file_ptr = fopen(logfile, "w");
	if (!log_file_ptr)
	{
		printf("Unable to open log file!\n");
		return 1;
	}

	fprintf(log_file_ptr,"Processing SRTM with input file: %s and output file: %s \n", infile, outfile);


	//Print time stamp
	time_t ltime;
	ltime = time(NULL);
	fprintf(log_file_ptr,"Starting process at: %s\n",asctime(localtime(&ltime)));


	//Allocate heap memory
	unsigned static long long block_pixels = IMAGEXDIM * (IMAGEYDIM / BLOCKS);
	unsigned static long long block_x_dim = IMAGEXDIM;
	unsigned static long long block_y_dim = IMAGEYDIM / BLOCKS;
	unsigned long long subblock_y_dim = block_y_dim / SUBBLOCKS;
	unsigned static long long out_block_x_dim = IMAGEXDIM / 3;
	unsigned static long long out_block_y_dim = (IMAGEYDIM / BLOCKS) / 3;
	unsigned long long out_block_pixels = out_block_x_dim * out_block_y_dim;
	unsigned long long subblock_pixels = block_pixels / SUBBLOCKS;
	unsigned long long out_subblock_pixels = subblock_pixels / 9;

	fprintf(log_file_ptr,"Allocating memory... input block is %llu pixels, output block is %llu pixels \n", block_pixels, out_block_pixels);
	int16_t* input_block = malloc(block_pixels * sizeof(int16_t));
	if (input_block == 0) {
		fprintf(log_file_ptr, "Unable to allocate memory for input block.  Exiting program.\n");
		fflush(log_file_ptr);
		return 1;
	}

	char* output_block = malloc(out_block_pixels);  //outblock is type byte, so sizeof returns 1
	if (output_block == 0) {
		fprintf(log_file_ptr, "Unable to allocate memory for output block.  Exiting program.\n");
		fflush(log_file_ptr);
		return 1;
	}

	//Open file for input
	fprintf(log_file_ptr,"Opening input file...\n");
	FILE* input_file_ptr;
	input_file_ptr = fopen(infile, "rb");

	if (!input_file_ptr)
	{
		fprintf(log_file_ptr,"Unable to open input file!\n");
		return 1;
	}

	//Open output file
	fprintf(log_file_ptr,"Opening output file...\n");
	fflush(log_file_ptr);
	FILE* output_file_ptr;
	output_file_ptr = fopen(outfile, "wb");
	if (!output_file_ptr)
	{
		fprintf(log_file_ptr,"Unable to open output file!\n");
		fflush(log_file_ptr);	
		return 1;
	}

	//Flush buffer
	fflush(stdout);
	fflush(log_file_ptr);

	//Processing loop
	static float vert_pixel_size = 30.87; //Vertical pixel size does not change.  30.87 meters
	static float top_latitude = 76.0; //Top latitude at the top of the entire image

	for(int block_i=0; block_i < BLOCKS; ++block_i) {
		ltime = time(NULL);
		fprintf(log_file_ptr,"Working on block %d of %d: %s \n", block_i+1, BLOCKS, asctime(localtime(&ltime)));
		fflush(log_file_ptr);	
		//Read input block
		fprintf(log_file_ptr,"Reading input block...\n");
		fflush(log_file_ptr);	
		fread(input_block, sizeof(int16_t), block_pixels, input_file_ptr);
		//Process this block, OpenMP parallelization happens here.
		fprintf(log_file_ptr,"Processing block...\n");
		fflush(log_file_ptr);	
		#pragma omp parallel
		{
			#pragma omp for
			for(unsigned long long subblock_i=0; subblock_i<SUBBLOCKS; ++subblock_i) {
				//Each subblock will start at the left most column, pointted to by subblock_ptr
				int16_t* subblock_ptr = input_block+(subblock_i*subblock_pixels);
				char* out_subblock_ptr = output_block+(subblock_i*out_subblock_pixels);
				float latitude, horz_pixel_size;
				float subblock_top_latitude = top_latitude - ((float)((block_i * block_y_dim)+subblock_i*subblock_y_dim))/60.0/60.0;   //Convert to degrees

				//Now, we can treat input and output subblock as standalone image (memory spacing works out)
				//As long as we use subblock_ptr and out_subblock_ptr and make sure we don't overrun 
				
				//Perform calculation of subblock
				for(unsigned long long j=0; j<out_block_y_dim/SUBBLOCKS; ++j) {
					//Calculate current latitude and horizontal pixel size at given j 
					latitude = subblock_top_latitude - (((float)j) + 0.5)/20.0/60.0;
					horz_pixel_size = 30.87 * cos(latitude*PI/180.0); // convert latitude to radians
					for(unsigned long long i=0; i<out_block_x_dim; ++i) {  //Subblocking is along y-axis only
						//Calculate maximum slope in degrees for the given output pixel
						int16_t* calc_block_ulpixel = subblock_ptr+(j*3ULL*IMAGEXDIM)+i*3ULL;
						*(out_subblock_ptr+j*out_block_x_dim+i) = maxdegreeslope(calc_block_ulpixel, 3, 3, IMAGEXDIM, horz_pixel_size, vert_pixel_size);

					}
				}
			} // End of subblock calculations
		}  //End of pragma omp parallel

		//At this point, the output block should all be calculated, write to file
		fprintf(log_file_ptr,"Writing output block...\n\n");
		fflush(log_file_ptr);	
		fwrite(output_block, sizeof(char), out_block_pixels, output_file_ptr);
	}//End of all blocks

	//Close files
	fclose(input_file_ptr);
	fclose(output_file_ptr);

	ltime = time(NULL);
	fprintf(log_file_ptr,"Done! %s\n", asctime(localtime(&ltime)));


}
