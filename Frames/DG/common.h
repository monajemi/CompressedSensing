#ifndef COMMON_H_
#define COMMON_H_

/* integer data types */
typedef unsigned int uint32_t;
typedef struct complex32 {
	int r;
	int i;
} complex32_t;

/* floating point data types */
typedef struct dblcomplex {
	double r;
	double i;
} dblcomplex_t;
dblcomplex_t dblcomplex_new(double r, double i);
double dblcomplex_norm_square(dblcomplex_t z);
dblcomplex_t dblcomplex_conjugate(dblcomplex_t z);
dblcomplex_t dblcomplex_plus(dblcomplex_t z1, dblcomplex_t z2);
dblcomplex_t dblcomplex_minus(dblcomplex_t z1, dblcomplex_t z2);
dblcomplex_t dblcomplex_times(dblcomplex_t z1, dblcomplex_t z2);
dblcomplex_t dblcomplex_zero();
dblcomplex_t dblcomplex_i_pow(uint32_t x);
void dblcomplex_print(dblcomplex_t z);

/* boolean data type */
typedef enum {
	FALSE = 0,
	TRUE = 1
} bool_t;

#define MAX_LINE_SIZE 30

/* Exponent of conjugate of i^k (0 -> 0, 1 -> 3, 2 -> 2, 3 -> 1, mod 4) */
#define I_CONJUGATE(k) ((((k) % 2) == 0) ? (k) : (((k) + 2) % 4))
#define PRINT_COMPLEX(z) if ((z).r == 0 && (z).i == 0) printf("0"); else if ((z).r == 0) printf("%di", (z).i); else if ((z).i == 0) printf("%d", (z).r); else printf("%d%s%di", (z).r, (z).i > 0 ? "+" : "", (z).i)
#define COMPLEX_NORM(z) ((z).r)

/* 2 ^ i */
#define POW(i) (1 << (i))
/* Is i odd or even? Assumption: working with unsigned integers */
#define ODD(i) ((i) & 1)
#define EVEN(i) !((i) & 1)

#define BITS_PER_BYTE 8
#define LEFT_MASK 0x80
#define RIGHT_MASK 0x01
#define ONE_MASK 0xff

#define ONES_MASK_32 0xffffffff
/* Get num_bits rightmost bits (all other bits set to zero) starting offset bits from right */
#define GET_BITS_32(val, num_bits, offset) (((val) >> (offset)) & (ONES_MASK_32 ^ (ONES_MASK_32 << (num_bits))))

/* Here arr is a char array, and we want to either get or set the ith bit in this
   array, where the zeroeth bit is the leftmost bit of the zeroeth character */
#define BIT_ARR_SET(arr, i, val) arr[(i) / BITS_PER_BYTE] = (arr[(i) / BITS_PER_BYTE] & (ONE_MASK ^ (0x01 << ((i) % BITS_PER_BYTE)))) | ((val) << ((i) % BITS_PER_BYTE))
#define BIT_ARR_GET(arr, i) ((arr[(i) / BITS_PER_BYTE] >> ((i) % BITS_PER_BYTE)) & 0x01)

/* Same as above, except each entry takes up two bytes (so can hold values 0,1,2,3) */
#define TWO_BIT_ARR_SET(arr, i, val) arr[(i) * 2 / BITS_PER_BYTE] = (arr[(i) * 2 / BITS_PER_BYTE] & (ONE_MASK ^ (0x03 << (((i) * 2) % BITS_PER_BYTE)))) | ((val) << (((i) * 2) % BITS_PER_BYTE))
#define TWO_BIT_ARR_GET(arr, i) ((arr[(i) * 2 / BITS_PER_BYTE] >> (((i) * 2) % BITS_PER_BYTE)) & 0x03)

#endif
