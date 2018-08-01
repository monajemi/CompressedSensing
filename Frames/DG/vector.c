#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "dgs.h"
/*#include "minpq.h"*/
#include "vector.h"

struct vector {
	uint32_t n;
	dblcomplex_t *vals;
};

vector_t vector_new(uint32_t n)
{
	vector_t v;
	
	v = (vector_t) malloc(sizeof(struct vector));
	v->n = n;
	v->vals = (dblcomplex_t *) malloc(sizeof(dblcomplex_t) * n);
	
	return v;
}

/* Create vector from image of k-sparse signal under image of d. */
vector_t vector_sample(dgs_t d, k_signal_t sig)
{
	uint32_t i, x;
	uint32_t N;
	vector_t v;
	dblcomplex_t *image_vals;
	dblcomplex_t weight, mat_entry;
	
	N = dgs_num_rows(d);
	v = vector_new(N);
	
	/* Add contributions of each relevant column from dg sensing matrix */
	
	image_vals = (dblcomplex_t *) malloc(sizeof(dblcomplex_t) * N);
	
	for (x = 0; x < N; x++)
		image_vals[x] = dblcomplex_zero();
	
	for (i = 0; i < sig->k; i++)
	{
		weight = sig->vals[i];
		
		for (x = 0; x < N; x++)
		{
			mat_entry = dblcomplex_i_pow(dgs_get_val(d, x, sig->indices[i]));
			image_vals[x] = dblcomplex_plus(image_vals[x], dblcomplex_times(weight, mat_entry));
		}
	}
	
	vector_set(v, image_vals);
	free(image_vals);
	
	return v;
}

vector_t vector_sample_transpose(dgs_t d, vector_t v)
{
	uint32_t i, j;
	uint32_t N, C;
	vector_t result;
	
	N = dgs_num_rows(d);
	C = dgs_num_cols(d);
	
	result = vector_new(C);
	
	for (j = 0; j < C; j++)
	{
		dblcomplex_t sum;
		
		sum = dblcomplex_zero();
		for (i = 0; i < N; i++)
			sum = dblcomplex_plus(dblcomplex_times(dblcomplex_conjugate(dblcomplex_i_pow(dgs_get_val(d, i, j))), v->vals[j]), sum);
		
		result->vals[j] = sum;
      
	}
	
	return result;
}

/* Not at all random right now */
vector_t vector_random_unit(dgs_t d, uint32_t k)
{
	uint32_t i;
	dblcomplex_t scratch;
	k_signal_t sig;
	
	sig = (k_signal_t) malloc(sizeof(struct k_signal));
	
	sig->k = k;
	sig->indices = (uint32_t *) malloc(sizeof(uint32_t) * k);
	sig->vals = (dblcomplex_t *) malloc(sizeof(dblcomplex_t) * k);
	
	for (i = 0; i < k; i++)
	{
		scratch.r = i + 2;
		scratch.i = 0;
		sig->vals[i] = dblcomplex_times(dblcomplex_i_pow(0), scratch);
		sig->indices[i] = (7 + i) * (i + 1);
	}
	
	return vector_sample(d, sig);
}

void vector_set(vector_t v, dblcomplex_t vals[])
{
	uint32_t i;
	
	for (i = 0; i < v->n; i++)
		v->vals[i] = vals[i];
	
	return;
}

void hadamard_helper(dblcomplex_t *vals, uint32_t m)
{
	uint32_t i;
	dblcomplex_t temp;

	if (m == 0)
		return;
	
	for (i = 0; i < POW(m - 1); i++)
	{
		temp = vals[i];
		vals[i] = dblcomplex_plus(vals[i], vals[i + POW(m - 1)]);
		vals[i + POW(m - 1)] = dblcomplex_minus(temp, vals[i + POW(m - 1)]);
	}
	
	hadamard_helper(vals, m - 1);
	hadamard_helper(&vals[POW(m - 1)], m - 1);
	
	return;
}

/* In place Hadamard transform of vector v */
void hadamard(vector_t v)
{
	uint32_t m, N;	/* N = 2^m = dim(v) */
	
	N = v->n;
	
	m = 0;
	while ((N / POW(m)) != 1)
		m++;
	
	hadamard_helper(v->vals, m);
	
	return;
}

/* Nodes and compare function for priority queue */
typedef struct node {
	double val;
	uint32_t index;
} *node_t;

node_t node_new(double val, uint32_t index)
{
	node_t n;
	
	n = (node_t) malloc(sizeof(struct node));
	
	n->val = val;
	n->index = index;
	
	return n;
}

void node_delete(node_t n)
{
	free(n);
}

/*
   compar(a,b) =  0 if a = b
			   =  1 if a > b
			   = -1 if a < b
*/

int compar_val(const void *a, const void *b)
{
	node_t *n1;
	node_t *n2;
	
	n1 = (node_t*) a;
	n2 = (node_t*) b;
	
	if ((*n1)->val > (*n2)->val)
		return 1;
	else if ((*n1)->val < (*n2)->val)
		return -1;
	
	return 0;
}

int compar_index(const void *a, const void *b)
{
	node_t n1;
	node_t n2;
	
	n1 = (node_t) a;
	n2 = (node_t) b;
	
	if (n1->index > n2->index)
		return 1;
	else if (n1->index < n2->index)
		return -1;
	
	return 0;
}

/* reconstruct k-sparse signal whose measurement through d is f */
int vector_reconstruct(vector_t f, uint32_t k, dgs_t d, k_signal_t sig)
{
    int C = dgs_num_cols(d);
   return  vector_large_reconstruct(f, k, d, sig, C);
	
}

int vector_large_reconstruct(vector_t f, uint32_t k, dgs_t d, k_signal_t sig, uint32_t Cm)
{
	uint32_t i, x, a, delta;
	uint32_t N, C;	/* dimensions of sensing matrix */
	vector_t f_shift;
	dblcomplex_t *lambda;
	dblcomplex_t *shift_vals;
	node_t *lambda_arr;
	
	N = dgs_num_rows(d);
	C = dgs_num_cols(d);
	
	/* Check that dimensions of sensing matrix agree with dimension of measurement v */
	if (N != f->n)
		return -1;
	
	/* Initialize array of peaks (4 values stored per peak) */
	lambda = (dblcomplex_t *) malloc(sizeof(dblcomplex_t) * Cm * 4);
	for (i = 0; i < Cm * 4; i++)
		lambda[i] = dblcomplex_zero();
	
	/* Calculate the peaks by averaging over offsets a in F_2^m */
	shift_vals = (dblcomplex_t *) malloc(sizeof(dblcomplex_t) * N);
	f_shift = vector_new(N);
	for (a = 0; a < N; a++)
	{	
		/* Calculate f(x xor a) * conj(f(x)) for all x */
		for (x = 0; x < N; x++)
			shift_vals[x] = dblcomplex_times(f->vals[x ^ a], dblcomplex_conjugate(f->vals[x]));
		vector_set(f_shift, shift_vals);
		
		/* Take Hadamard transform of resulting vector (in place) */
		hadamard(f_shift);
		
		/* For each peak delta, calculate lambda(delta, a), and add to lamda(delta) */
		for (delta = 0; delta < Cm; delta++)
		{
			uint32_t lam_index;	/* based off i^(-a(Q_delta)aT) */
			dblcomplex_t fourier_coef;	/* entry of hadamard(f_shift) indexed by vector a(Q_delta) */
			
			lam_index = 4 * delta + dgs_get_val(d, a, delta);
			fourier_coef = f_shift->vals[dgs_get_vec(d, a, delta)];
			
			lambda[lam_index] = dblcomplex_plus(lambda[lam_index], fourier_coef);
		}
	}
	vector_delete(f_shift);
	free(shift_vals);
	
	/* Find k highest peaks (support of signal) */
	lambda_arr = (node_t*) malloc(sizeof(node_t) * Cm);
	
	for (i = 0; i < Cm; i++)
	{
		dblcomplex_t val;
		node_t n;
		
		val = dblcomplex_zero();
		val = dblcomplex_plus(val, dblcomplex_times(lambda[i * 4], dblcomplex_new(1, 0)));
		val = dblcomplex_plus(val, dblcomplex_times(lambda[i * 4 + 1], dblcomplex_new(0, -1)));
		val = dblcomplex_plus(val, dblcomplex_times(lambda[i * 4 + 2], dblcomplex_new(-1, 0)));
		val = dblcomplex_plus(val, dblcomplex_times(lambda[i * 4 + 2], dblcomplex_new(0, 1)));
		
		n = node_new(dblcomplex_norm_square(val), i);
		
		lambda_arr[i] = n;
	}
	free(lambda);
	
	
	qsort(lambda_arr, Cm, sizeof(node_t), &compar_val);
	
	for (i = 0; i < k; i++)
		sig->indices[i] = lambda_arr[Cm - i - 1]->index;
	
	for (i = 0; i < Cm; i++)
		node_delete(lambda_arr[i]);
	free(lambda_arr);
	
	/* TODO: Retrieve k signal entries */
	
	return 0;
}

void vector_print(vector_t v)
{
	uint32_t i;
	
	for (i = 0; i < v->n; i++)
	{
		dblcomplex_print(v->vals[i]);
		printf("\n");
	}
	printf("\n");
}

void vector_delete(vector_t v)
{
	free(v->vals);
	free(v);
}
