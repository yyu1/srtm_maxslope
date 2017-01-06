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
#define IMAGEXDIM 1296000
#define IMAGEYDIM 475200
#define BLOCKS 10
#define SUBBLOCKS 20



int main() {
	
	char infile[100] = "/nobackupp6/nexprojects/CMS-ALOS/srtm1sec_mosaic/global_srtm3_aster_mosaic_1sec.int"
	char outfile[100] = "/nobackupp6/nexprojects/CMS-ALOS/srtm1sec_mosaic/global_srtm3_aster_maxdegreeslope_3sec.byt"

	printf("Processing SRTM with input file: %s and output file: %s \n", infile, outfile);


	//Print time stamp
	time_t ltime;
	ltime = time(NULL);
	printf("Starting process at: %s\n",asctime(localtime(&ltime)));


	//Allocate heap memory
	unsigned static long long block_pixels = IMAGEXDIM * (IMAGEYDIM / BLOCKS);
	unsigned static long long block_x_dim = IMAGEXDIM;
	unsigned static long long block_y_dim = IMAGEYDIM / BLOCKS;
	unsigned static long long subblock_y_dim = block_y_dim / SUBBLOCKS;
	unsigned static long out_block_x_dim = IMAGEXDIM / 3;
	unsigned static long out_block_y_dim = (IMAGEYDIM / BLOCKS) / 3;
	unsigned static long long out_block_pixels = out_block_x_dim * out_block_y_dim;

	printf("Allocating memory...\n");
	int16_t* input_block = malloc(block_pixels * sizeof(int16_t));
	char* output_block = malloc(out_block_pixels);  //outblock is type byte, so sizeof returns 1

	//Open file for input
	printf("Opening input file...\n");
	FILE* input_file_ptr;
	input_file_ptr = fopen(infile, "rb");

	if (!input_file_ptr)
	{
		printf("Unable to open input file!\n");
		return 1;
	}

	//Open output file
	printf("Opening output file...\n");
	FILE* output_file_ptr;
	output_file_ptr = fopen(outfile, "wb");
	if (!output_file_ptr)
	{
		printf("Unable to open output file!\n");
		return 1;
	}

	//Processing loop
	static float vert_pixel_size = 30.87 //Vertical pixel size does not change.  30.87 meters
	static float top_latitude = 76.0 //Top latitude at the top of the entire image

	for(int block_i=0; block_i < BLOCKS; ++block_i) {
		ltime = time(NULL);
		printf("Working on block %d of %d: %s \n", block_i+1, BLOCKS, asctime(localtime(&ltime)));
		
		//Read input block
		printf("Reading input block...\n");
		fread(input_block, sizeof(int16_t), block_pixels, input_file_ptr);

		//Process this block, OpenMP parallelization happens here.
		printf("Processing block...\n");
		#pragma omp parallel
		{
			#pragma omp for
			for(int subblock_i=0; subblock_i<SUBBLOCKS; ++subblock_i) {
				//Each subblock will start at the left most column, pointted to by subblock_ptr
				int16_t* subblock_ptr = input_block+(subblock_i*subblock_pixels);
				char* out_subblock_ptr = output_block+(subblock_i*subblock_pixels);
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
						*(out_subblock_ptr+j*out_block_xdim+i) = maxdegreeslope(subblock_ptr, 3, 3, IMAGEXDIM, horz_pixel_size, vert_pixel_size);

					}
				}
			} // End of subblock calculations
		}  //End of pragma omp parallel

		//At this point, the output block should all be calculated, write to file
		printf("Writing output block...\n\n");
		fwrite(output_block, sizeof(char), out_block_pixels, output_file_ptr);
	}//End of all blocks

	//Close files
	fclose(input_file_ptr);
	fclose(output_file_ptr);

	ltime = time(NULL);
	printf("Done! %s\n", asctime(localtime(&ltime)));


}
