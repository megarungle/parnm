// Copyright Utkin Konstantin, 2022

#define DEBUG 1

#include <math.h>
#include <memory.h>

#if DEBUG
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#endif

void PrintMatrix(double* M, int n, int m);
void PullMatrix(double* A, double* B, int n, int m, int x1, int y1, int x2, int y2);
void Cholesky_Decomposition(double* A, double* L, int n);
void _Cholesky_Decomposition(double* A, double* L, int n);
void TriangularSystemSolve(double* L21, double* L11, double* A21, int n, int m);

#if DEBUG
int main()
{
    int n = 7;
    int r = 3;

    double* A = new double[n * n];
    
    srand(time(NULL));
    
    // Generate test matrix
    for (int i = 0; i < n; ++i)
    {
        for (int j = i; j < n; ++j)
        {
            A[i * n + j] = double(rand() % 100) + 1;
            A[j * n + i] = A[i * n + j];
        }
    }

    // A > 0
    for (int i = 0; i < n; ++i)
    {
        double sum = 0;
        A[i * n + i] = 0;
        for (int j = 0; j < n; ++j)
        {
            sum += A[i * n + j];
        }
        A[i * n + i] = sum;
    }

    PrintMatrix(A, n, n);

    double* A11 = new double[r * r];
    double* A12 = new double[(n - r) * r];
    double* A21 = new double[r * (n - r)];
    double* A22 = new double[(n - r) * (n - r)];

    PullMatrix(A, A11, n, r,     0, 0, r, r);
    PullMatrix(A, A12, n, n - r, 0, r, r, n);
    PullMatrix(A, A21, n, r,     r, 0, n, r);
    PullMatrix(A, A22, n, n - r, r, r, n, n);

    /*
    PrintMatrix(A11, r, r);
    PrintMatrix(A12, n - r, r);
    PrintMatrix(A21, r, n - r);
    PrintMatrix(A22, n - r, n - r);
    */

    double* L11 = new double[r * r];
    _Cholesky_Decomposition(A11, L11, r);

    PrintMatrix(L11, r, r);

    double* L21 = new double[r * (n - r)];
    TriangularSystemSolve(L21, L11, A21, r, n - r);

    PrintMatrix(L21, r, n - r);

    delete[] A;
    //delete[] L;

    return 0;
}
#endif

#if DEBUG
void PrintMatrix(double* M, int n, int m)
{
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            printf("%6.3f ", M[i * n + j]);
        }
        printf("\n");
    }
    printf("\n");
}
#endif

void PullMatrix(double* A, double* B, int n, int m, int x1, int y1, int x2, int y2)
{
    for (int i = 0; i < x2 - x1; ++i)
    {
        for (int j = 0; j < y2 - y1; ++j)
        {
            B[i * m + j] = A[x1 * n + y1 + i * n + j];
        }
    }
}

void Cholesky_Decomposition(double* A, double* L, int n)
{
    int r = 10;

    // L11
    _Cholesky_Decomposition(A, L, r);

    // L21

    // L22 recursion
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

void TriangularSystemSolve(double* L21, double* L11, double* A21, int n, int m)
{
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double sum = 0.0;
            for (int k = 0; k < n; ++k)
            {
                if (j == 0)
                {
                    L21[i * n + j] = A21[i * n + j] / L11[0];
                    break;
                }
                else
                {
                    if (k != j)
                    {
                        sum += L21[i * n + k] * L11[j * n + k];  // L11 transp
                    }
                    else
                    {
                        L21[i * n + k] = (A21[i * n + k] - sum) / L11[j * n + k];  // L11 transp
                        sum = 0.0;
                        break;
                    }
                }
            }
        }
    }
}
