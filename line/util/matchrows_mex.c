#include "mex.h"
void mexFunction(int nl, mxArray *l[], int nr, const mxArray *r[])
{
    double *A, *B;
    A  = mxGetPr(r[0]);
    B  = mxGetPr(r[1]);
    mwSize n, m, i, j, k, t, c;
    n = mxGetM(r[0]); // rows A
    m = mxGetN(r[0]); // columns A
    c = mxGetM(r[1]); // rows b
    double *v;
    v = (double *) mxMalloc( c * sizeof(*v) );
    for (t = 0; t < c; t++) {
        v[t] = (double) -1.0;
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                if (A[i+j*n] != B[t+c*j]) {
                    break;
                }
            }
            if (j == m) {
                v[t] = ((double) (i + 1));
                break;
            }
        }
    }
    l[0] = mxCreateDoubleMatrix(c,1,mxREAL);
    mxSetPr(l[0], v);
}
