/* File:       conj_grad.c
 * Author:     Vincent Zhang
 *
 * Purpose:    A serial conjugate gradient solver program. Due to time limits,
 *             the MPI parallel version is not included in this source code.
 *
 * Compile:    gcc -g -Wall -lm -o conj_grad conj_grad.c
 * Run:        conj_grad [order] [tolerance] [iterations]
 *                       [Optional suppress output(n)] < [file]
 *
 * Input:      A file that contains a symmetric, positive definite matrix A,
 *             and the corresponding right hand side vector B. Preferably, each
 *             line consists of [n] elements and the [n+1] line would be the b.
 *
 * Output:     1. The number of iterations,
 *             2. The time used by the solver (not including I/O),
 *             3. The solution to the linear system (if not suppressed),
 *             4. The norm of the residual calculated by the conjugate gradient
 *                method, and
 *             5. The norm of the residual calculated directly from the
 *                definition of residual.
 *				6. if debug compiled a file name 4compute  will be created. It contains
 *				the system to be solved. It can be used to check the result with wolframalpha.

 * Algorithm:  The matrix A's initially read and parsed into an one-dimensional
 *             array; the right hand side vector b is stored in an array as
 *             well. After some preparation work of allocating memory and
 *             assigning variables the program jumps into the main loop, the
 *             conjugate gradient solver. For the exact mathematical procedure,
 *             please refer to http://www.cs.usfca.edu/~peter/cs625/prog2.pdf
 *             and http://en.wikipedia.org/wiki/Conjugate_gradient_method for a
 *             much better demonstration.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// #include <mpi.h>     /* For MPI functions, etc */
#include "timer.h"
#include "generate_equation.h"

/**
 * Calculates the dot product of two arrays.
 *
 * @param a The first array.
 * @param b The second array.
 * @param size The size of the arrays.
 * @return The dot product of the two arrays.
 */
double dotProduct(double *a, double *b, int size)
{
	double sum = 0.0;
	int i;
	for (i = 0; i < size; i++)
	{
		sum += (a[i] * b[i]);
		// printf("Dot product check, round: %d, sum = %lf\n", i, sum);
	}
	// printf("==== Dot product check, final: %lf\n", sum);
	return sum;
}

/**
 * Multiplies each element of a vector by a scalar value.
 *
 * @param dest The destination vector where the result will be stored.
 * @param v The input vector to be multiplied.
 * @param s The scalar value to multiply each element of the vector by.
 * @param size The size of the vector.
 * @return The destination vector with the scalar-vector product.
 */
double *scalarVector(double *dest, double *v, double s, int size)
{
	// printf("== Begin: scalar vector product ==");
	int i;
	for (i = 0; i < size; i++)
	{
		dest[i] = s * v[i];
		// printf("Scalar vector product check, round: %d, sum = %lf\n", i, dest[i]);
	}
	return dest;
}

/**
 * Adds two vectors element-wise and stores the result in the destination vector.
 *
 * @param dest The destination vector where the result will be stored.
 * @param a The first input vector.
 * @param b The second input vector.
 * @param size The size of the vectors.
 * @return The destination vector with the element-wise sum of the input vectors.
 */
double *vectorAdd(double *dest, double *a, double *b, int size)
{
	// printf("== Begin: vector adding vector ==");
	int i;
	for (i = 0; i < size; i++)
	{
		dest[i] = a[i] + b[i];
		// printf("Vector add vector check, round: %d, sum = %lf\n", i, dest[i]);
	}
	return dest;
}

/**
 * Subtracts two vectors element-wise and stores the result in the destination vector.
 *
 * @param dest The destination vector to store the result.
 * @param a The first vector.
 * @param b The second vector.
 * @param size The size of the vectors.
 * @return The destination vector with the subtracted elements.
 */
double *vectorSubtract(double *dest, double *a, double *b, int size)
{
	// printf("== Begin: vector subtracting vector ==");
	int i;
	for (i = 0; i < size; i++)
	{
		dest[i] = a[i] - b[i];
		// printf("Vector subtract vector check, round: %d, sum = %lf\n", i, dest[i]);
	}
	return dest;
}

/**
 * Performs matrix-vector multiplication.
 *
 * This function multiplies a matrix by a vector and stores the result in the destination array.
 *
 * @param dest   The destination array to store the result.
 * @param matrix The matrix to be multiplied.
 * @param v      The vector to be multiplied.
 * @param size   The size of the matrix and vector.
 * @return       The destination array containing the result of the multiplication.
 */
double *matrixVector(double *dest, double *matrix, double *v, int size)
{
	/* The Cross Product */
	// printf("== Begin: matrix vector product ==");
	int i, j;
	for (i = 0; i < size; i++)
	{
		dest[i] = 0.0;
		for (j = 0; j < size; j++)
		{
			dest[i] += matrix[i * size + j] * v[j];
			// printf("Matrix vector check, round: %d, sum = %lf\n", i, dest[i]);
		}
	}
	return dest;
}

/**
 * @brief Assigns the elements of one array to another array.
 *
 * This function copies the elements from array 'b' to array 'a' of size 'size'.
 *
 * @param a Pointer to the destination array.
 * @param b Pointer to the source array.
 * @param size The number of elements to copy.
 */
void assignVector(double *a, double *b, int size)
{
	/* Essentially does the same job as memcpy. */
	int i;
	for (i = 0; i < size; i++)
	{
		a[i] = b[i];
	}
}

/**
 * @file conj_grad.c
 * @brief Conjugate Gradient Solver
 *
 * This program solves a linear system of equations using the Conjugate Gradient method.
 * It takes command line arguments for the dimension of the matrix, tolerance, maximum iterations,
 * and an optional argument to suppress output.
 *
 * The program generates a random linear system, reads the matrix and right-hand side vector from file,
 * and performs the Conjugate Gradient iterations until the tolerance is reached or the maximum iterations
 * are exceeded. The solution vector is then printed.
 *
 * @param argc The number of command line arguments
 * @param argv An array of strings containing the command line arguments
 * @return 0 if the program executed successfully, 1 otherwise
 */
int main(int argc, char **argv)
{
	int i, /*j,*/ k, order, max_iterations, suppress_output;
	double start, finish, elapsed;
	double tolerance;
	generateRandomLinearSystem(1);
	// for (i = 0; i < argc; i++) {
	// 	  printf("##### %d = %s\n", i, argv[i]);
	// }

	suppress_output = 0;

	if (argc > 3)
	{
		/* Parse command line arguments */
		order = atoi(argv[1]); /* Determines dimension of matrix A */
							   // printf("==== Dimension of matrix: %d\n", order);

		tolerance = atof(argv[2]);

		max_iterations = atoi(argv[3]);

		/* Use string compare to check optional suppress */
		if (argc == 5 && (strcmp(argv[4], "n") == 0))
		{
			printf("==== %s\n", argv[4]);
			suppress_output = 1;
		}
		else
		{
			printf("==== No optional command \n");
		}

		/* 1-d array for A */
		double *matrix = generateRandomLinearSystem(order);
		/* Right hand side vector, B */
		double *rhs = malloc(order * sizeof(double));

		/* Read matrix (A) from file and put into array of doubles */
		// for (i = 0; i < order; i++)
		// { /* number of rows */
		// 	for (j = 0; j < order; j++)
		// 	{ /* number of cols */
		// 		if (!scanf("%lf", &matrix[(i * order) + j]))
		// 		{
		// 			break;
		// 		}
		// 	}
		// }

		/* Read right hand side (B) */
		rhs = solutionVector(order);

		printf("==== PRINTING MATRIX ==== \n");
		int newline = 0;
		for (i = 0; i < (order * order); i++)
		{
			if (newline == order)
			{
				newline = 0;
				printf("\n");
			}
			printf("%f ", matrix[i]);
			newline++;
		}
		printf("\n==== PRINTING MATRIX ==== \n");

		printf("==== PRINTING RIGHT HAND SIDE ==== \n");
		for (i = 0; i < order; i++)
		{
			printf("%f\n", rhs[i]);
		}
		printf("==== PRINTING RIGHT HAND SIDE ==== \n");

#ifdef DEBUG1
		/*print matrix on file*/
		FILE *file = fopen("4compute.txt", "w");
		if (file == NULL)
		{
			printf("Impossibile aprire il file.\n");
			return 1;
		}

		for (int i = 0; i < order; i++)
		{
			for (int j = 0; j < order; j++)
			{
				fprintf(file, "%.3lfx_%d", matrix[i * order + j], j + 1);
				if (j < order - 1)
				{
					fprintf(file, "+");
				}
				else
				{
					fprintf(file, "= %.3lf", rhs[i]);
					if (i < order - 1)
					{
						fprintf(file, ",");
					}
				}
			}
		}

		fclose(file);

		/*end print*/
#endif
		/* Prepare variables for the main loop */
		k = 0;
		double beta;
		double alpha;
		double *x = malloc(order * sizeof(double));
		double *s = malloc(order * sizeof(double));
		double *x_prev = malloc(order * sizeof(double));
		double *p = malloc(order * sizeof(double));
		double *p_prev = malloc(order * sizeof(double));
		double *r = malloc(order * sizeof(double));
		double *r_prev = malloc(order * sizeof(double));
		double *r_prev_prev = malloc(order * sizeof(double));
		double *holderVector = malloc(order * sizeof(double));

		/* Before we start, copy values of B into R_0 */
		// memcpy(r, rhs, (order * sizeof(double)));
		// memcpy(r_prev, rhs, (order * sizeof(double)));
		// memcpy(r_prev_prev, rhs, (order * sizeof(double)));
		for (i = 0; i < order; i++)
		{
			x[i] = 0;
			r[i] = rhs[i];
			r_prev[i] = rhs[i];
			r_prev_prev[i] = rhs[i];
		}

		GET_TIME(start);

		while ((k < max_iterations) && (dotProduct(r, r, order) > tolerance))
		{
			// memcpy(r_prev_prev, r_prev, (order * sizeof(double)));
			assignVector(r_prev_prev, r_prev, order);
			assignVector(r_prev, r, order);
			assignVector(p_prev, p, order);
			assignVector(x_prev, x, order);
			k++;
			if (k == 1)
			{

				/* P_1 = R_0 */
				memcpy(p, r_prev, (order * sizeof(double)));
				memcpy(p_prev, p, (order * sizeof(double)));
				// assignVector(p, r_prev, order);
				// assignVector(p_prev, p, order);
			}
			else
			{
				/* BETA_k = [R_(k-1) * R_(k-1)] / [R_(k-2) * R_(k-2)]  */
				beta = dotProduct(r_prev, r_prev, order) / dotProduct(r_prev_prev, r_prev_prev, order);

				// memcpy(p_prev, p, (order * sizeof(double)));

				/* P_k = R_(k-1) + [BETA_k * P_(k-1)] */
				holderVector = scalarVector(holderVector, p_prev, beta, order);
				p = vectorAdd(p, r, holderVector, order);
			}
			/* S_k = (A * P_k) */
			s = matrixVector(s, matrix, p, order);

			// memcpy(r_prev, r, (order * sizeof(double)));

			/* ALPHA_k = [R_(k-1) * R_(k-1)] / [P_k * S_k] */
			double d1 = dotProduct(r_prev, r_prev, order);
			double d2 = dotProduct(p, s, order);
			alpha = d1 / d2;

			// memcpy(x_prev, x, (order * sizeof(double)));

			/* X_k = X_(k-1) + (ALPHA_k * P_k) */
			holderVector = scalarVector(holderVector, p, alpha, order);
			x = vectorAdd(x, x_prev, holderVector, order);

			/* R_k = R_(k-1) - (ALPHA_k * S_k) */
			holderVector = scalarVector(holderVector, s, alpha, order);
			r = vectorSubtract(r, r_prev, holderVector, order);
		}

		GET_TIME(finish);
		elapsed = finish - start;

		printf("========= Solver Completed ========= \n");
		printf("The code to be timed took %lf seconds\n", elapsed);
		printf("Number of iterations: %d \n", k);

		if (suppress_output == 0)
		{
			printf("Solution to the matrix: \n");
			for (i = 0; i < (order); i++)
			{
				printf("%f\n", x[i]);
			}
		}

		/* Free all objects that allocated memory once at the end */
		free(matrix);
		free(rhs);

		free(x);
		free(s);
		free(p);
		free(r);

		free(x_prev);
		free(p_prev);
		free(r_prev);
		free(r_prev_prev);
		free(holderVector);
	}
	else
	{
		printf("Usage: %s [order] [tolerance] [iterations] [Optional suppress output(n)]\n", argv[0]);
	}

	return 0;
} /* main a*/
