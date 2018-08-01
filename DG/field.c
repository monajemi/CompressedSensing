#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "common.h"
#include "field.h"

struct field {
	uint32_t deg;		/* F_2^deg */
	char *poly;			/* Generating primitive polynomial */
	char *trace_by_exp;	/* table of trace values for nonzero field elements, [e] = Tr(a^e) */
	char *trace_by_vec; /* table of trace values for all field elements (x_0, ..., x_m-1) */
};

/* helper function for field_new */
void field_generate(field_t f);

field_t field_new(char *s)
{
	uint32_t i;
	uint32_t deg;
	field_t f;
	
	deg = 0;
	for (i = 0; i < strlen(s) && isdigit(s[i]); i++)
		deg = deg * 10 + (s[i] - '0');
	
	if (i == strlen(s) || s[i] != '\t')
		return NULL;
	
	if (EVEN(deg) || deg < 3)
		return NULL;
	
	f = (field_t) malloc(sizeof(struct field));
	f->deg = deg;
	f->poly = (char*) malloc(deg + 1);
	
	s += (i + 1);
	for (i = 0; i < strlen(s) && i <= deg; i++)
		f->poly[i] = s[i] - '0';
	
	if (i <= deg)
	{
		free(f->poly);
		free(f);
		return NULL;
	}
	
	field_generate(f);
	
	return f;
}

void field_delete(field_t f)
{
	free(f->poly);
	free(f->trace_by_exp);
	free(f);
}

uint32_t p_to_uint(uint32_t len, char *p)
{
	uint32_t i;
	uint32_t result;
	
	result = 0;
	for (i = 0; i < len; i++)
		result += POW(i) * p[i];
	
	return result;
}

int p_right_shift(uint32_t len, char *p)
{
	uint32_t i;
	int result;
	
	result = p[len - 1];
	
	for (i = len - 1; i > 0; i--)
		p[i] = p[i - 1];
	
	p[0] = 0;
	
	return result;
}

void p_add(uint32_t len, char *p1, char *p2)
{
	uint32_t i;
	
	for (i = 0; i < len; i++)
		p2[i] ^= p1[i];
}

char *p_clone(uint32_t len, char *p)
{
	uint32_t i;
	char *pclone;
	
	pclone = (char*) malloc(len);
	
	for (i = 0; i < len; i++)
		pclone[i] = p[i];
	
	return pclone;
}

void field_generate(field_t f)
{
	uint32_t e, i;
	char **rep_table;
	
	f->trace_by_exp = (char*) malloc(POW(f->deg) - 1);
	
	/* Generate m-tuple representation of all field elements a^e */
	rep_table = (char**) malloc(sizeof(char*) * (POW(f->deg) - 1));
	rep_table[0] = (char*) calloc(f->deg, 1);
	rep_table[0][0] = 1;
	for (e = 1; e < POW(f->deg) - 1; e++)
	{
		rep_table[e] = p_clone(f->deg, rep_table[e - 1]);
		if (p_right_shift(f->deg, rep_table[e]) == 1)
			p_add(f->deg, f->poly, rep_table[e]);
	}
	
	/* Find trace of a^0, a^1, ..., a^(m-1) */
	for (e = 0; e < f->deg; e++)
	{
		char *p;
		
		p = (char*) calloc(f->deg, 1);
		
		for (i = 0; i < f->deg; i++)
			p_add(f->deg, rep_table[(e * POW(i)) % (POW(f->deg) - 1)], p);
		
		f->trace_by_exp[e] = p[0];
		
		free(p);
	}
	
	/* Find trace of a^m, ..., a^(2^m - 2) using trace values of standard basis */
	for (e = f->deg; e < POW(f->deg) - 1; e++)
	{
		f->trace_by_exp[e] = 0;
		for (i = 0; i < f->deg; i++)
		{
			if (rep_table[e][i] == 1)
				f->trace_by_exp[e] ^= f->trace_by_exp[i];
		}
	}
	
	/* Fill in trace by m-tuple representation */
	f->trace_by_vec = (char*) malloc(POW(f->deg));
	f->trace_by_vec[0] = 0;
	for (e = 0; e < POW(f->deg) - 1; e++)
		f->trace_by_vec[p_to_uint(f->deg, rep_table[e])] = f->trace_by_exp[e];
	
	for (e = 0; e < POW(f->deg) - 1; e++)
		free(rep_table[e]);
	free(rep_table);
}

uint32_t field_get_degree(field_t f)
{
	return f->deg;
}

/* Return the trace of a^e */
char field_get_trace(field_t f, uint32_t e)
{
	return f->trace_by_exp[e % (POW(f->deg) - 1)];
}

char field_get_trace_by_vec(field_t f, uint32_t v)
{
	return f->trace_by_vec[v];
}

void field_print(field_t f)
{
	uint32_t i;
	
	/* Print trace values */
	printf("%d: 0", f->deg);
	for (i = 0; i < POW(f->deg) - 1; i++)
		printf("%d", (int) field_get_trace(f, i));
	printf("\n");
}

int field_verify(field_t f)
{
	uint32_t e;
	uint32_t i, j;
	uint32_t deg;
	bool_t non_degenerate;
	deg = field_get_degree(f);
	
	/* Ensure non-degeneracy of trace (at least one element with nonzero trace) */
	
	non_degenerate = FALSE;
	for (e = 0; e < POW(deg) - 1; e++)
	{
		if (field_get_trace(f, e) == 1)
			non_degenerate = TRUE;
	}
	
	if (!non_degenerate)
	{
		printf("Trace zero for all field elements: not non-degenerate\n");
		return -1;
	}
	/*mexPrintf("Sina");*/
	
	/* Ensure linearity of trace */
	for (i = 0; i < POW(deg); i++)
	{
		for (j = i; j < POW(deg); j++)
		{
			if (field_get_trace_by_vec(f, i ^ j) != (field_get_trace_by_vec(f, i) ^ field_get_trace_by_vec(f, j)))
			{
				printf("Trace nonlinear!\n");
				printf("Tr(%d) = %d\n", i, field_get_trace_by_vec(f, i));
				printf("Tr(%d) = %d\n", j, field_get_trace_by_vec(f, j));
				printf("Tr(%d) = %d\n", i ^ j, field_get_trace_by_vec(f, i ^ j));
				return -1;
			}
		}
	}
	
	return 0;
}
