#ifndef FIELD_H_
#define FIELD_H_

typedef struct field *field_t;

field_t field_new(char *s);
void field_delete(field_t f);

uint32_t field_get_degree(field_t f);
char field_get_trace(field_t f, uint32_t e);
char field_get_trace_by_vec(field_t f, uint32_t v);

int field_verify(field_t f);
void field_print(field_t f);

#endif
