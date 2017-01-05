#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.1415929653589

//This function returns the value at the pixel given by i,j in the block with the given pointer to the upper left pixel
inline int16_t value(int16_t* ulpixel, int i, int j, long imagehorzsize) {
	return *(ulpixel + (j*imagehorzsize) + i);
}


//ulpixel is a pointer to the upper left pixel of the block
//blockhorzsize is the horizontal dimension of the block to be calculated
//blockvertsize is the vertical dimension of the block to be calculated
//imagehorzsize is the horizontal dimension of the overall input image, used to calculate pointers for pixel in the block
//horzdistance is the distance between horizontal pixels, in the same unit as the DEM
//vertdistance is the distance between vertical pixels, in the same unit as the DEM
char maxdegreeslope(int16_t* ulpixel, int blockhorzsize, int blockvertsize, long imagehorzsize, float horzdistance, float vertdistance)  {

	float maxslope = 0;
	float tempslope;

	int secondi, secondj;  //i and j of the other location for calculating slope

	for (int i=0; i<blockhorzsize; ++i) {
		for (int j=0; j<blockvertsize; ++j) {

			//This part of the loop iterates over each pixel inside the block

			//left
			secondi = i-1;
			secondj = j;
			if (secondi >= 0) {
				tempslope = fabsf((float)value(ulpixel,i,j,imagehorzsize) - (float)value(ulpixel,secondi,secondj,imagehorzsize))/horzdistance;
				if (tempslope > maxslope) { maxslope = tempslope; }
			}

			//right
			secondi = i+1;
			secondj = j;
			if (secondi < blockhorzsize) {
				tempslope = fabsf((float)value(ulpixel,i,j,imagehorzsize) - (float)value(ulpixel,secondi,secondj,imagehorzsize))/horzdistance;
				if (tempslope > maxslope) { maxslope = tempslope; }
			}

			//lower left
			secondi = i-1;
			secondj = j+1;
			if (secondi >=0 && secondj < blockvertsize) {
				tempslope = fabsf((float)value(ulpixel,i,j,imagehorzsize) - (float)value(ulpixel,secondi,secondj,imagehorzsize))/sqrt(horzdistance*horzdistance + vertdistance * vertdistance);
				if (tempslope > maxslope) { maxslope = tempslope; }
			}

			//lower
			secondi = i;
			secondj = j+1;
			if (secondj < blockvertsize) {
				tempslope = fabsf((float)value(ulpixel,i,j,imagehorzsize) - (float)value(ulpixel,secondi,secondj,imagehorzsize))/vertdistance;
				if (tempslope > maxslope) { maxslope = tempslope; }
			}


			//lower right
			secondi = i+1;
			secondj = j+1;
			if (secondi < blockhorzsize && secondj < blockvertsize) {
				tempslope = fabsf((float)value(ulpixel,i,j,imagehorzsize) - (float)value(ulpixel,secondi,secondj,imagehorzsize))/sqrt(horzdistance*horzdistance + vertdistance * vertdistance);
				if (tempslope > maxslope) { maxslope = tempslope; }
			}

		}
	}

	return atan(maxslope)*180/PI;  //chnage to degree slope
}
