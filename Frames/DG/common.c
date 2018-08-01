#include <stdio.h>
#include "common.h"

dblcomplex_t dblcomplex_new(double r, double i)
{
	dblcomplex_t z;
	
	z.r = r;
	z.i = i;
	
	return z;
}

double dblcomplex_norm_square(dblcomplex_t z)
{
	return z.r * z.r + z.i * z.i;
}

dblcomplex_t dblcomplex_plus(dblcomplex_t z1, dblcomplex_t z2)
{
	dblcomplex_t result;
	
	result.r = z1.r + z2.r;
	result.i = z1.i + z2.i;
	
	return result;
}

dblcomplex_t dblcomplex_minus(dblcomplex_t z1, dblcomplex_t z2)
{
	dblcomplex_t result;
	
	result.r = z1.r - z2.r;
	result.i = z1.i - z2.i;
	
	return result;
}

dblcomplex_t dblcomplex_times(dblcomplex_t z1, dblcomplex_t z2)
{
	dblcomplex_t result;
	
	result.r = z1.r * z2.r - z1.i * z2.i;
	result.i = z1.r * z2.i + z1.i * z2.r;
	
	return result;
}

dblcomplex_t dblcomplex_conjugate(dblcomplex_t z)
{
	dblcomplex_t result;
	
	result.r = z.r;
	result.i = -z.i;
	
	return result;
}

dblcomplex_t dblcomplex_zero()
{
	dblcomplex_t result;
	
	result.r = 0.0;
	result.i = 0.0;
	
	return result;
}

dblcomplex_t dblcomplex_i_pow(uint32_t x)
{
	dblcomplex_t result;
	
	switch (x % 4)
	{
		case 0:
			result.r = 1.0;
			result.i = 0.0;
			break;
		case 1:
			result.r = 0.0;
			result.i = 1.0;
			break;
		case 2:
			result.r = -1.0;
			result.i = 0.0;
			break;
		case 3:
			result.r = 0.0;
			result.i = -1.0;
			break;
		default:
			printf("PROBLEM with dblcomplex_i_pow\n");
			break;
	}
	
	return result;
}

void dblcomplex_print(dblcomplex_t z)
{
	if ((z.r == 0.0) && (z.i == 0.0))
		printf("%.3f", z.r);
	else if (z.r == 0.0)
		printf("%.3fi", z.i);
	else if (z.i == 0.0)
		printf("%.3f", z.r);
	else
		printf("%.3f%s%.3fi ", z.r, (z.i < 0.0) ? "  " : " + ", z.i);
}




