/* Delsarte-Goethals sensing matrix */

#ifndef DGS_H_
#define DGS_H_

typedef struct dgs *dgs_t;

#include "field.h"
dgs_t dgs_new(field_t f, uint32_t r);
void dgs_delete(dgs_t d);

/* Basic information about matrix */
uint32_t dgs_get_m(dgs_t d);
uint32_t dgs_get_r(dgs_t d);
uint32_t dgs_num_rows(dgs_t d);
uint32_t dgs_num_cols(dgs_t d);
uint32_t dgs_get_val(dgs_t d, uint32_t x, uint32_t P);	/* x is row index, P column index; returns xPxT */
uint32_t dgs_get_vec(dgs_t d, uint32_t x, uint32_t P);	/* returns binary vector xP */

/* Return information about inner products of columns */
complex32_t dgs_inner_prod(dgs_t d, uint32_t col_1, uint32_t col_2);
int dgs_inner_prod_norm(dgs_t d, uint32_t col_1, uint32_t col_2);

/* Print stuff about matrix */
void dgs_print(dgs_t d);					/* Print matrix */
void dgs_print_col(dgs_t d, uint32_t col);
void dgs_print_inner_prod(dgs_t d); 		/* Print matrix formed from inner products of column vectors */
void dgs_print_inner_prod_norm(dgs_t d);	/* Print matrix formed from norms of inner products of column vectors */

int dgs_verify(dgs_t d);	/* unit test sensing matrix */

#endif
