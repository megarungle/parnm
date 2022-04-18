#define DEBUG 1

#include <math.h>
#include <memory.h>

#if DEBUG
#include <stdio.h>
#endif

#if DEBUG
void main()
{
    
}
#endif

void Cholesky_Decomposition(double* A, double* L, int n)
{

}

void _Cholesky_Decomposition(double* A, double* L, int n)
{
    memset(L, 0, sizeof(double) * n * n);

    for (int i = 0; i < n; ++i)
    {
        L[i * n + i] = A[i * n + i];
        for (int k = 0; k < i; ++k)
        {
            L[i * n + i] -= L[i * n + k] * L[i * n + k];
        }

        L[i * n + i] = sqrt(L[i * n + i]);

        for (int j = i + 1; j < n; ++j)
        {
            L[j * n + i] = A[j * n + i];
            for (int k = 0; k < i; ++k)
            {
                L[j * n + i] -= L[i * n + k] * L[j * n + k];
            }

            L[j * n + i] /= L[i * n + i];
        }
    }
}
