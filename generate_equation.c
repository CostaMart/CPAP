#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// #include <mpi.h>     /* For MPI functions, etc */
#include <time.h>
// verifica sia simmetrica
#include "generate_equation.h"
#include <stdbool.h>

/**
 * Checks if a given matrix is positive symmetric.
 *
 * @param matrix The matrix to be checked.
 * @param n The size of the matrix.
 * @return True if the matrix is positive symmetric, false otherwise.
 */
bool isPositiveSymmetric(double **matrix, int n)
{
    // Verifica la simmetria
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (matrix[i][j] != matrix[j][i])
            {
                return false;
            }
        }
    }
    return true;
}

/**
 * Generates a medium matrix by averaging the elements of two input matrices and adding the order to the diagonal.
 *
 * @param m1 The first input matrix.
 * @param m2 The second input matrix.
 * @param order The order of the matrices.
 * @return The resulting medium matrix.
 */
double **mediumMatrix(double **m1, double **m2, int order)
{
    double **result = malloc(order * sizeof(double *));

    for (int i = 0; i < order; i++)
    {
        result[i] = malloc(order * sizeof(double));
        for (int j = 0; j < order; j++)
        {
            result[i][j] = 0.5 * (m1[i][j] + m2[i][j]);
        }
    }

    // aggiungo NI alla diagonale
    for (int i = 0; i < order; i++)
    {
        result[i][i] += order;
    }

#ifdef DEBUG
    printf("\n La matrice è simmetrica positiva? %s\n", isPositiveSymmetric(result, order) ? "Sì" : "No");
#endif
    return result;
}

/**
 * Generates a random linear system of equations.
 *
 * @param order The order of the linear system.
 * @return A pointer to the generated linear system.
 */
double *generateRandomLinearSystem(int order)
{
    srand(time(NULL));
    double **A = malloc(order * sizeof(double *));
    double **AT = malloc(order * sizeof(double *));
    double random;
    double *toreturn;

    // genera le righe per la matrice
    for (int i = 0; i < order; i++)
    {
        A[i] = malloc(order * sizeof(double));
        AT[i] = malloc(order * sizeof(double));
    }

    // popola la matrice di random double
    for (int i = 0; i < order; i++)
    {
        for (int j = 0; j < order; j++)
        {
            random = (double)rand() / RAND_MAX;
            A[i][j] = random;
            AT[j][i] = random;
        }
    }

#ifdef DEBUG
    for (int i = 0; i < order; i++)
    {
        printf("\n");
        for (int j = 0; j < order; j++)
        {
            printf("%f ", A[i][j]);
        }
    }

    printf("\n");
    printf("\n");
    for (int i = 0; i < order; i++)
    {
        printf("\n");
        for (int j = 0; j < order; j++)
        {
            printf("%f ", AT[i][j]);
        }
    }
#endif

    double **result;
    result = mediumMatrix(A, AT, order);

#ifdef DEBUG
    printf("\n");
    printf("\n");
    for (int i = 0; i < order; i++)
    {
        printf("\n");
        for (int j = 0; j < order; j++)
        {
            printf("%f ", result[i][j]);
        }
    }
#endif

    // converti la matrice ottenuta in una rappresentazione 1D perchè il programma principale utilizza questa
    toreturn = malloc(order * order * sizeof(double));
    int n = 0;
    for (int i = 0; i < order; i++)
    {
        for (int j = 0; j < order; j++)
        {
            toreturn[n] = result[i][j];
            n++;
        }
    }

#ifdef DEBUG
    printf("\n");
    for (int i = 0; i < order * order; i++)
    {
        printf("|| %f ||", toreturn[i]);
    }
#endif
    free(AT);
    free(A);
    return toreturn;
}

/**
 * Generates a random solution vector.
 *
 * @param length The length of the solution vector.
 * @return The generated solution vector.
 */
double *solutionVector(int length)
{
    // Generate the solution vector B
    double *B = malloc(length * sizeof(double));
    for (int i = 0; i < length; i++)
    {
        B[i] = (double)rand() / RAND_MAX;
    }

#ifdef DEBUG
    printf("\n vettore soluzioni: \n[");
    for (int i = 0; i < length; i++)
    {
        printf("%f, ", B[i]);
    }
    printf("]\n");
#endif
    return B;
}
