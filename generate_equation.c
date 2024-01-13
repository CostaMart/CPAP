#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// #include <mpi.h>     /* For MPI functions, etc */
#include <time.h>
// verifica sia simmetrica
#include "generate_equation.h"
#include <stdbool.h>

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

    // TODO: Verifica che tutti gli autovalori siano positivi.
    // Questo è un problema più complesso che richiede la decomposizione
    // in autovalori di una matrice, che va oltre il codice di base in C.

    return true;
}

// ---------- Calcola la matrice media
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

//-------- genera un sistema lienare con matrice dei coefficienti simmetrica, così da poter applicare il gradiente coniugato
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

// ------ genera un vettore di soluzioni randomico
double *solutionVector(int lunghezza)
{
    // abbiamo la matrice ora generiamo un b
    double *B = malloc(lunghezza * sizeof(double));
    for (int i = 0; i < lunghezza; i++)
    {
        B[i] = (double)rand() / RAND_MAX;
    }

    printf("\n vettore soluzioni: \n[");
    for (int i = 0; i < lunghezza; i++)
    {
        printf("%f, ", B[i]);
    }
    printf("]\n");
    return B;
}
