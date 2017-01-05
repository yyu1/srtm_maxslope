//Testing maxdegreeslope function
#include "maxdegreeslope.h"
#include <stdio.h>

int main() {
	int16_t block[9] = {5,10,12,6,10,8,7,3,9};	

	printf("max slope is: %d\n", maxdegreeslope(&(block[0]), 3, 3, 3, 10.0, 10.0));

	return 0;
}
