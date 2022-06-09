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
void CholeskySeq(double* A, double* L, int n);
void TriangularSystemSolve(double* L21, double* L11, double* A21, int n, int m);
void _Cholesky_Decomposition(double* A, double* L, int n, int r, int iter, int N);
void CholeskyFactor(double* A22, double* L21, double* A_22, int n, int r);

#if DEBUG
int main()
{
    int n = 300;

    double* A = new double[n * n];
    double* L = new double[n * n];

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

    Cholesky_Decomposition(A, L, n);

    PrintMatrix(L, n, n);

    return 0;
}
#else
int main()
{
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
    memset(L, 0, sizeof(double) * n * n);

    int r = n;

    if (n > 20)
    {
        r = 20;
        if (n > 50)
        {
            r = 50;
            if (n > 100)
            {
                r = 100;
            }
        }
    }

    _Cholesky_Decomposition(A, L, n, r, 0, n);
}

void _Cholesky_Decomposition(double* A, double* L, int n, int r, int iter, int N)
{
    if (n == r)
    {
        double* A22 = A;

        double* L22 = new double[n * n];
        CholeskySeq(A22, L22, n);

#if DEBUG
        printf("Iter: %d, L22:\n", iter);
        PrintMatrix(L22, n, n);
#endif

        // Push to L
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                L[n * iter + N * n * iter + i * N + j] = L22[i * n + j];
            }
        }

        delete[] L22;
        delete[] A22;

        return;
    }

    // A11 working
    double* A11 = new double[r * r];
    PullMatrix(A, A11, n, r, 0, 0, r, r);

    double* L11 = new double[r * r];
    CholeskySeq(A11, L11, r);

#if DEBUG
    printf("Iter: %d, L11:\n", iter);
    PrintMatrix(L11, r, r);
#endif


    // A21 working
    double* A21 = new double[r * (n - r)];
    PullMatrix(A, A21, n, r, r, 0, n, r);

    double* L21 = new double[r * (n - r)];
    TriangularSystemSolve(L21, L11, A21, r, n - r);
    
#if DEBUG
    printf("Iter: %d, L21:\n", iter);
    PrintMatrix(L21, r, n - r);
#endif


    // A22 recursive working
    if (n > r)
    {
        double* A22 = new double[(n - r) * (n - r)];
        PullMatrix(A, A22, n, n - r, r, r, n, n);

        double* A_22 = new double[(n - r) * (n - r)];
        CholeskyFactor(A22, L21, A_22, n - r, r);

        delete[] A22;

        _Cholesky_Decomposition(A_22, L, n - r, r, iter + 1, N);
    }
    else
    {
        int size = n - r;
        int offset = r;

        double* A22 = new double[size * size];
        PullMatrix(A, A22, n, n - r, r, r, n, n);

        double* L22 = new double[size * size];
        CholeskySeq(A22, L22, size);

#if DEBUG
        printf("Iter: %d, L22:\n", iter);
        PrintMatrix(L22, size, size);
#endif

        // Push to L
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                L[offset * iter + N * offset * iter + i * N + j] = L22[i * size + j];
            }
        }

        delete[] L22;
        delete[] A22;
    }


    // L building
    // L11
    for (int i = 0; i < r; ++i)
    {
        for (int j = 0; j < r; ++j)
        {
            L[r * iter + N * r * iter + i * N + j] = L11[i * r + j];
        }
    }
    // L21
    for (int i = 0; i < n - r; ++i)
    {
        for (int j = 0; j < r; ++j)
        {
            L[r * iter + N * r * iter + (i + r) * N + j] = L21[i * r + j];
        }
    }

    delete[] A;
    delete[] A11;
    delete[] A21;
    delete[] L11;
    delete[] L21;
}

void CholeskySeq(double* A, double* L, int n)
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
    memset(L21, 0, sizeof(double) * n * m);

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

void CholeskyFactor(double* A22, double* L21, double* A_22, int n, int r)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double mult = 0.0;

            for (int k = 0; k < r; ++k)
            {
                mult += L21[i * r + k] * L21[j * r + k];
            }

            A_22[i * n + j] = A22[i * n + j] - mult;

            mult = 0.0;
        }
    }
}
