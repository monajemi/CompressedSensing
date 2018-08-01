#ifndef VECTOR_H
#define VECTOR_H

typedef struct vector *vector_t;

typedef struct k_signal {
	uint32_t k;
	uint32_t *indices;
	dblcomplex_t *vals;
} *k_signal_t;

vector_t vector_new(uint32_t n);
vector_t vector_sample(dgs_t d, k_signal_t sig);
vector_t vector_sample_transpose(dgs_t d, vector_t v);
vector_t vector_random_unit(dgs_t d, uint32_t k);
void vector_set(vector_t v, dblcomplex_t vals[]);

int vector_large_reconstruct(vector_t v, uint32_t k, dgs_t d, k_signal_t sig,uint32_t C_m);
int vector_reconstruct(vector_t v, uint32_t k, dgs_t d, k_signal_t sig);

void vector_print(vector_t v);

void vector_delete(vector_t v);

#endif
