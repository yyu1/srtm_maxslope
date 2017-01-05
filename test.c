//Testing maxdegreeslope function
#include "maxdegreeslope.h"
#include <stdio.h>

int main() {
	int16_t block[16] = {0,0,0,0,0,5,1,12,0,6,10,8,0,7,3,9};	

	printf("max slope is: %d\n", maxdegreeslope(&(block[8]), 2, 2, 4, 10.0, 10.0));

	return 0;
}
