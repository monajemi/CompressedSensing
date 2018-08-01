#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "field.h"
#include "dgs.h"

struct dgs {
	uint32_t m;
	uint32_t r;
	char ***p_val;
	short ***p_vec;
};

/* helper functions for constructing DG sensing matrix */
char **calc_dg_p_val(field_t f, uint32_t t);
short **calc_dg_p_vec(field_t f, uint32_t t);

dgs_t dgs_new(field_t f, uint32_t r)
{
	dgs_t d;
	uint32_t m;
	uint32_t t;
	
	if (f == NULL)
		return NULL;
	
	m = field_get_degree(f);
	
	/* m should be odd and >= 3 */
	if ((m % 2 == 0) || (m < 3))
		return NULL;
	
	/* r should be between 0 and (m-1)/2 */
	if (r > (m - 1) / 2)
		return NULL;
	
	d = (dgs_t) malloc(sizeof(struct dgs));
	
	d->m = m;
	d->r = r;
	
	d->p_val = (char ***) malloc(sizeof(char**) * (r + 1));
	for (t = 0; t <= r; t++)
	{
		d->p_val[t] = calc_dg_p_val(f, t);
	}
	d->p_vec = (short ***) malloc(sizeof(short**) * (r + 1));
	for (t = 0; t <= r; t++)
	{	
        d->p_vec[t] = calc_dg_p_vec(f, t);
	}
	return d;
}

void dgs_delete(dgs_t d) {

	uint32_t i, j;
	
	for (i = 0; i <= d->r; i++)
	{
		for (j = 0; j < d->m; j++)
		{
			free(d->p_val[i][j]);
			free(d->p_vec[i][j]);
		}
		free(d->p_val[i]);
		free(d->p_vec[i]);
	}
	
	free(d->p_val);
	free(d->p_vec);

	free(d);
}

/* Get entry (i, j) for Kerdock matrix P0(a^e) (this function is not used for zero matrix) */
uint32_t kerdock_entry(field_t f, uint32_t e, uint32_t i, uint32_t j)
{
	return field_get_trace(f, e + i + j);
}

/* Get entry (i, j) for matrix Pt(a^e) (function not used for zero matrix) */
uint32_t pt_entry(field_t f, uint32_t t, uint32_t e, uint32_t i, uint32_t j)
{
	uint32_t val1;
	uint32_t val2;
	
	val1 = field_get_trace(f, i + j * POW(t) + e);
	val2 = field_get_trace(f, i * POW(t) + j + e);
	
	return val1 ^ val2;
}

/* IMPORTANT: x-ordering must be shared between calc_kerdock_p_val and calc_dg_p_val.
   Also, all values calculated for following helper functions are mod 4. */

/* Calculate xPxT for all matrices Pt(a^e), for all x */
char **calc_dg_p_val(field_t f, uint32_t t)
{
	uint32_t e, x, i, j;
	uint32_t degree;
	char **p_val;
	
	degree = field_get_degree(f);
	p_val = (char**) malloc(sizeof(char*) * POW(degree));
	
	/* Entry x for xP0(a^e)xT */
	
	for (x = 0; x < POW(degree); x++)
	{
		p_val[x] = (char*) malloc(POW(degree));
		
		/* Zero matrix results in zero for all values of x */
		p_val[x][0] = 0;
		
		for (e = 0; e < POW(degree) - 1; e++)
		{
			uint32_t prod;
			
			if (t == 0)
			{
				/* Calculate xP0(a^e)xT */
				prod = 0;
				for (i = 0; i < degree; i++)
				{
					char x_i = (x >> (degree - i - 1)) & RIGHT_MASK;
					
					if (x_i != 0)
					{
						prod += kerdock_entry(f, e, i, i) * x_i;
						for (j = i + 1; j < degree; j++)
						{
							char x_j = (x >> (degree - j - 1)) & RIGHT_MASK;		
							prod += 2 * x_i * x_j * kerdock_entry(f, e, i, j);
						}
					}
				}
			}
			else
			{
				/* Calculate xPt(a^e)xT for t > 0*/
				prod = 0;
				for (i = 0; i < degree; i++)
				{
					char x_i = (x >> (degree - i - 1)) & RIGHT_MASK;
					
					if (x_i != 0)
					{
						/* No diagonal! */
						for (j = i + 1; j < degree; j++)
						{
							char x_j = (x >> (degree - j - 1)) & RIGHT_MASK;
							prod += 2 * x_i * x_j * pt_entry(f, t, e, i, j);
						}
					}
				}
			}
			
			p_val[x][e + 1] = prod % 4;
		}
	}
	
	return p_val;
}

/* Calculate xP (mod 2), stored as a binary vector, for all matrices Pt(a^e), for all x */
short **calc_dg_p_vec(field_t f, uint32_t t)
{
	uint32_t e, x, i, j;
	uint32_t degree;
	uint32_t char_arr_size;
	short **p_vec;
	
	degree = field_get_degree(f);
	p_vec = (short**) malloc(sizeof(short*) * POW(degree));
	
	/* Entry x for xPt(a^e)xT ranging over e */
	
	char_arr_size = POW(degree) * sizeof(short);
	
	for (x = 0; x < POW(degree); x++)
	{
		p_vec[x] = (short*) malloc(char_arr_size);
		
		/* Zero matrix */
		
		p_vec[x][0] = 0;
		
		for (e = 0; e < POW(degree) - 1; e++)
		{
			uint32_t vec;
			
			/* Calculate xPt(a^e) */
			vec = 0;
			for (j = 0; j < degree; j++)
			{
				uint32_t prod;
				
				if (t == 0)
				{
					/* prod = (xP0(a^e))_j (mod 2) */
					prod = 0;
					for (i = 0; i < degree; i++)
						prod += ((x >> (degree - i - 1)) & RIGHT_MASK) * kerdock_entry(f, e, i, j);
					prod %= 2;
				}
				else
				{
					/* prod = (xPt(a^e))_j (mod 2) for t > 0 */
					prod = 0;
					for (i = 0; i < degree; i++)
						prod += ((x >> (degree - i - 1)) & RIGHT_MASK) * pt_entry(f, t, e, i, j);
					prod %= 2;
				}
				
				vec |= prod * POW(degree - j - 1);
			}
			
			p_vec[x][e + 1] = vec;
		}
	}
	
	return p_vec;
}

complex32_t dgs_inner_prod(dgs_t d, uint32_t col_1, uint32_t col_2)
{
	uint32_t x, t;
	uint32_t *p_indices_1;
	uint32_t *p_indices_2;
	complex32_t result;
	
	/* Rightmost m bits define P0, next m bits define P1, next m bits define P2, ... */
	
	p_indices_1 = (uint32_t *) malloc(sizeof(uint32_t) * (d->r + 1));
	p_indices_2 = (uint32_t *) malloc(sizeof(uint32_t) * (d->r + 1));
	
	for (t = 0; t <= d->r; t++)
	{
		p_indices_1[t] = GET_BITS_32(col_1, d->m, t * d->m);
		p_indices_2[t] = GET_BITS_32(col_2, d->m, t * d->m);
	}
		
	result.r = 0;
	result.i = 0;
	
	for (x = 0; x < POW(d->m); x++)
	{
		uint32_t p1, p2;
		
		p1 = 0;
		p2 = 0;
		for (t = 0; t <= d->r; t++)
		{
			p1 += TWO_BIT_ARR_GET(d->p_val[t][p_indices_1[t]], x);
			p2 += TWO_BIT_ARR_GET(d->p_val[t][p_indices_2[t]], x);
		}
		
		/* Recall that complex inner product uses conjugates of values in second vector */
		p2 = I_CONJUGATE(p2);
		
		switch((p1 + p2) % 4)
		{
			case 0:
				result.r++;
				break;
			case 1:
				result.i++;
				break;
			case 2:
				result.r--;
				break;
			case 3:
				result.i--;
				break;
			default:
				printf("ERRORORRORORO!!!\n");
				break;
		}
	}
	
	return result;
}

int dgs_inner_prod_norm(dgs_t d, uint32_t col_1, uint32_t col_2)
{
	complex32_t z;
	
	z = dgs_inner_prod(d, col_1, col_2);
	
	return z.r * z.r + z.i * z.i;
}

void dgs_print_col(dgs_t d, uint32_t col)
{
	uint32_t x;
	
	for (x = 0; x < POW(d->m); x++)
	{	
		switch (dgs_get_val(d, x, col))
		{
			case 0:
				printf(" 1  ");
				break;
			case 1:
				printf(" i  ");
				break;
			case 2:
				printf("-1  ");
				break;
			case 3:
				printf("-i  ");
				break;
			default:
				printf("BIG PROBLEM\n");
				break;
		}
	}
	
	printf("\t");
	
	for (x = 0; x < POW(d->m); x++)
		printf("%d ", dgs_get_vec(d, x, col));
	
	printf("\n");
}

void dgs_print(dgs_t d)
{
	uint32_t p;
	
	/* Range over all matrices in DG(m, r) */
	for (p = 0; p < POW(d->m * (d->r + 1)); p++)
		dgs_print_col(d, p);
}

void dgs_print_inner_prod(dgs_t d)
{
	uint32_t i, j;
	complex32_t inner_prod;

	for (i = 0; i < POW(d->m * (1 + d->r)); i++)
	{
		for (j = 0; j < POW(d->m * (1 + d->r)); j++)
		{
			inner_prod = dgs_inner_prod(d, i, j);
			PRINT_COMPLEX(inner_prod);
			printf(" ");
		}
		printf("\n");
	}
}

void dgs_print_inner_prod_norm(dgs_t d)
{
	uint32_t i, j;
	
	for (i = 0; i < POW(d->m * (1 + d->r)); i++)
	{
		for (j = 0; j < POW(d->m * (1 + d->r)); j++)
			printf("%d ", dgs_inner_prod_norm(d, i, j));
		printf("\n");
	}
}

uint32_t dgs_num_rows(dgs_t d)
{
	return POW(d->m);
}

uint32_t dgs_num_cols(dgs_t d)
{
	return POW(d->m * (d->r + 1));
}

uint32_t dgs_get_val(dgs_t d, uint32_t x, uint32_t P)
{
	uint32_t result;
	uint32_t t;
	
	result = 0;
	for (t = 0; t <= d->r; t++)
		result += d->p_val[t][x][GET_BITS_32(P, d->m, d->m * t)];
	result %= 4;
	
	return result;
}

uint32_t dgs_get_vec(dgs_t d, uint32_t x, uint32_t P)
{
	uint32_t result;
	uint32_t t;
	
	result = 0;
	for (t = 0; t <= d->r; t++)
		result ^= d->p_vec[t][x][GET_BITS_32(P, d->m, d->m * t)];
	
	return result;
}

uint32_t dgs_get_m(dgs_t d)
{
	return d->m;
}

uint32_t dgs_get_r(dgs_t d)
{
	return d->r;
}


bool_t row_sums_zero(dgs_t d)
{
	uint32_t x, P, N, C;
	
	N = dgs_num_rows(d);
	C = dgs_num_cols(d);
	
	for (x = 0; x < N; x++)
	{
		int r, i;
		
		r = 0;
		i = 0;
		
		for (P = 0; P < C; P++)
		{
			switch (dgs_get_val(d, x, P))
			{
				case 0:
					r++;
					break;
				case 1:
					i++;
					break;
				case 2:
					r--;
					break;
				case 3:
					i--;
					break;
				default:
					break;
			}
		}
		
		if (r != 0 || i != 0)
		{
			printf("Not all row sums vanish!\n");
			printf("Row %d sums to %d + %di\n", x, r, i);
			return FALSE;
		}
	}
	
	return TRUE;
}

bool_t rows_orthogonal(dgs_t d)
{
	uint32_t x1, x2, P, N, C;
	
	N = dgs_num_rows(d);
	C = dgs_num_cols(d);
	
	for (x1 = 0; x1 < N; x1++)
	{
		for (x2 = x1 + 1; x2 < N; x2++)
		{
			int r, i;
			
			r = 0;
			i = 0;
			
			for (P = 0; P < C; P++)
			{
				switch((dgs_get_val(d, x1, P) + I_CONJUGATE(dgs_get_val(d, x2, P))) % 4)
				{
					case 0:
						r++;
						break;
					case 1:
						i++;
						break;
					case 2:
						r--;
						break;
					case 3:
						i--;
						break;
					default:
						break;
				}
			}
			
			if (r != 0 || i != 0)
			{
				printf("Not all rows are orthogonal!\n");
				printf("Rows %d and %d have inner product %d + %di\n", x1, x2, r, i);
				return FALSE;
			}
		}
	}
	
	return TRUE;
}

int dgs_verify(dgs_t d)
{
	uint32_t m, r, N, C;
	
	m = dgs_get_m(d);
	r = dgs_get_r(d);
	N = dgs_num_rows(d);
	C = dgs_num_cols(d);

	/* Check for consistency of m, r */
	
	if (EVEN(m))
	{
		printf("m = %d not odd!\n", m);
		return -1;
	}
	if (m < 3)
	{
		printf("m = %d too small!\n", m);
		return -1;
	}
	if (r > (m - 1) / 2)
	{
		printf("r = %d > (m - 1) / 2 = (%d - 1) / 2\n", r, m);
		return -1;
	}
	
	/* Check for consistency of N, C*/
	
	if (N != POW(m))
	{
		printf("N = %d != 2^m = 2^%d\n", N, m);
		return -1;
	}
	if (C != POW((r + 1) * m))
	{
		printf("C = %d != 2^((r + 1) * m) = 2^((%d + 1) * %d)\n", C, r, m);
		return -1;
	}
	
	/* Ensure rows are orthogonal */
	
	/*if (!rows_orthogonal(d))
		return -1;
	*/
	/* Ensure row sums are zero */
	/*
	if (!row_sums_zero(d))
		return -1;
	*/
	/* Ensure columns form group under pointwise multiplication */
	
	/* Ensure column sums are bounded */
	
	return 0;
}















