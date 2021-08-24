#include "mex.h"
void mexFunction(int nl, mxArray *l[], int nr, const mxArray *r[])
{
    double *A, *b;
    A  = mxGetPr(r[0]);
    b  = mxGetPr(r[1]);  
    mwSize n, m, i, j, k;
    n = mxGetM(r[0]);
    m = mxGetN(r[0]);
    for (i = 0; i < n; i++) {
       k = i;
       for (j = 0; j < m; j++) {
          if (A[k] != b[j]) {
             break;
          }
          k += n;
       }
       if (j == m) {
          l[0] = mxCreateDoubleScalar((double) (i + 1));
          return;
       }
    }
    l[0] = mxCreateDoubleScalar((double) -1);
}
