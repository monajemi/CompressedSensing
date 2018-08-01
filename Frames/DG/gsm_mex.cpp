/********************************************************************** 
 * gsm.cpp
 *generates the sensing matrix
 * phi=Phi(m,r,normalized)
 / *
 ********************************************************************/
#include <matrix.h>
#include <math.h>
#include <mex.h>   
#include "field.c"
#include "field.h"
#include "dgs.c"
#include "dgs.h"
#include "vector.h"
#include "vector.c"
#include "common.c"
#include "common.h"
//#include "test.h"
//#include "test.c"
//#include "main.c"

#define MAX_DEGREE 17

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//declare variables
    mxArray *m_in_m, *r_in_m, *normalized_in_m, *phi_out_m;
    const mwSize *dims;
    double  *phi, *normalized;
    double  *m,*r;
    int dimx, dimy, numdims,k ;
 
//associate inputs
     m_in_m = mxDuplicateArray(prhs[0]);
    r_in_m = mxDuplicateArray(prhs[1]);
    normalized_in_m = mxDuplicateArray(prhs[2]);

//figure out dimensions
//    dims = mxGetDimensions(prhs[0]);
  //  k=(int) dims[1];
   // numdims = mxGetNumberOfDimensions(prhs[0]);
    //dimy = (int)dims[0]; dimx = (int)dims[1];

//associate pointers
    m = mxGetPr(m_in_m);
    int N=(int)pow(2,*m);
    r = mxGetPr(r_in_m);
    int C=(int) pow(2,(*r+1)*(*m));
    normalized=mxGetPr(normalized_in_m);

//associate outputs
    phi_out_m = plhs[0] = mxCreateDoubleMatrix(N,C,mxREAL);    
    phi = mxGetPr(phi_out_m);
  
//generating the finite field:
   
   char s[MAX_LINE_SIZE];

	FILE *poly_stream;
	FILE *incremental;
	field_t *fields;
	dgs_t d1;
	
	fields = (field_t *) malloc(sizeof(field_t) );
	
	poly_stream = fopen("polynomials.txt", "r");

 		int mDG=(int) (*m);
		switch(mDG)
		{ 
			case 3:
				 strcpy(s,"3\t110100");
				break;
			case 5:
				strcpy(s,"5\t101001");
				break;
			case 7:
                strcpy(s,"7\t100100010");
                break;
			case 9:
                strcpy(s,"9\t100010000100");
                break;
			case 11:
                strcpy(s,"11\t101000000001");
                break;
			case 13:
                strcpy(s,"13\t110110000000010");
                break;
			case 15:
				strcpy(s,"15\t1100000000000001");
				break;
			default:
				strcpy(s,"5\t101001");
				break;
		}	
		field_t f;
		uint32_t degree;
		f = field_new(s);
		if (f != NULL)
		{
			degree = field_get_degree(f);
			
			if (ODD(degree) && (degree == *m) )
     				fields[0] = f;
           		 else
				field_delete(f);
		}

	d1 = dgs_new(fields[0], (int)*r);
     for(int j=0;j<C;j++)
        for(int i=0;i<N;i++)
            phi[j*N+i]=dgs_get_val(d1, i, j);
	
    field_delete(fields[0]);
    dgs_delete(d1);
   return;
}
