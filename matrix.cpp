#include "matrix.h"

void matrix::set_value(int row, int col, double value)
{
	int position = (row - 1) * (matrix::ncols) + col - 1;
	matrix::body[position] = value;
}

double matrix::get_value(int row, int col)
{
	int position = (row - 1) * (matrix::ncols)+col - 1;
	double ans = body[position];
	return ans;
}

void matrix::del_self()
{
	delete[] matrix::body;
}

matrix addition(matrix& mat1, matrix& mat2)
{
	if (mat1.nrows == mat2.nrows and mat1.ncols == mat2.ncols)
	{
		int size = mat1.ncols * mat1.nrows;
		matrix ans(mat1.nrows, mat1.ncols);
		for (size_t i = 0; i < size; i++)
		{
			ans.body[i] = mat1.body[i] + mat2.body[i];
		}
		return ans;
	}
	matrix ans(1, 1);
	return ans;
}

matrix matrix_multiplication(matrix& mat1, matrix& mat2)
{
	double val=0;
	double val1=0;
	double val2=0;

	if (mat1.ncols == mat2.nrows)
	{
		matrix ans(mat1.nrows, mat2.ncols);
		int size = mat1.nrows * mat2.nrows;

		for (int j = 0; j < mat1.nrows; j++)
		{
			for (int k = 0; k < mat2.ncols; k++)
			{
				val = 0;
				for (int i = 0; i < mat1.ncols; i++)
				{
					val1 = mat1.get_value(j+1, i+1);
					val2 = mat2.get_value(i+1, k+1);
					val = val + (val1 * val2);
				}
				ans.set_value(j+1, k+1, val);
			}
		}
		return ans;
	}
	else
	{
		return matrix(1, 1);
	}
	
}

matrix scalar_multiplication(matrix& mat, double scalar)
{
	matrix ans(mat.nrows, mat.ncols);
	int size = mat.nrows * mat.ncols;
	for (size_t i = 0; i < size; i++)
	{
		ans.body[i] = mat.body[i] * scalar;
	}
	return ans;
}

matrix transpose(matrix& mat)
{
	matrix ans(mat.ncols, mat.nrows);
	double val;
	for (int i = 0; i < mat.nrows; i++)
	{
		for (int j = 0; j < mat.ncols; j++)
		{
			val = mat.get_value(i+1, j+1);
			ans.set_value(j+1, i+1, val);
		}
	}
	return ans;
}

matrix solve(matrix& mat1, matrix& mat2)
{
	//mat1 should have size MxN. mat2 should have size Mx1. Solution will be size Mx1.
	
	return matrix(1,1);
}

int max_in_column(matrix& mat, int column, int start)
{
	double val1;
	double val2=0;
	int ans = 0;
	int n = mat.nrows + 1;

	for (size_t i = start; i < n; i++)
	{
		val1 = mat.get_value(i, column);
		val1 = abs(val1);
			if (val1 > val2)
			{
				val2 = val1;
				ans = i;
			};
	}
	return ans;
}

matrix make_pivot(matrix& mat, int row1, int row2) {
	//row1 and row2 must both be 1<= x <= mat.nrows

	// mat must be square
	if (mat.nrows != mat.ncols)
	{
		return mat;
	}

	else
	{
		for (int i = 1; i < mat.nrows + 1; i++)
		{
			for (int j = 1; j < mat.nrows + 1; j++)
			{
				if (i == j)
				{
					mat.set_value(i, j, 1);
				}
				else
				{
					mat.set_value(i, j, 0);
				}
			}
		}

		//if any row1 or row2 are <1 or >number of rows in mat or is row1 and row2 are equal then return an identity matrix the same size as mat.
		if (row1 < 1 or row1>mat.nrows or row2<1 or row2>mat.nrows or row1 == row2) {
			return mat;
		}

		//if all conditions are met then create a pivot matrix (Identity matrix with row 1 swapped with row 2)
		else
		{


			mat.set_value(row1, row1, 0);
			mat.set_value(row1, row2, 1);

			mat.set_value(row2, row2, 0);
			mat.set_value(row2, row1, 1);

			return mat;
		}
	}
}

matrix invert_pivot(matrix& mat) {
	matrix temp(mat.nrows, mat.ncols);
	for (int i = 1; i < temp.nrows + 1; i++)
	{
		for (int j = 1; j < temp.ncols + 1; j++)
		{
			temp.set_value(temp.nrows + 1 - i, temp.ncols + 1 - j, mat.get_value(i,j));
		}
	}

	return temp;
}

matrix forward_substitution(matrix& lower_trianglular, matrix& b) {
	matrix y(b.nrows, 1);
	double inner_sum;
	y.set_value(1, 1, b.get_value(1, 1));

	for (int i = 2; i < y.nrows + 1; i++)
	{
		inner_sum = 0;
		for (int j = 0; j < i; j++)
		{
			inner_sum = inner_sum + lower_trianglular.get_value(i, j)* y.get_value(j, 1);
		}
		y.set_value(i, 1, b.get_value(i, 1)-inner_sum);
	}

	return y;
}

matrix LU(matrix& mat) {
	//mat should be a square matrix (size MxM)
	int pivot;
	int n = mat.nrows;
	int N = n + 1;
	double scaling_factor = 1;
	double val = 0;
	
	matrix P(n, n);
	matrix L(n, n);
	matrix U(n, n);
	matrix M(n, n);
	matrix A(n, n);
	matrix c(n, 1);
	

	A = make_pivot(A, 0, 0);

	for (int i = 0; i < n*n; i++)
	{
		U.body[i] = mat.body[i];
	}

	L = make_pivot(L, 0, 0);
	for (int i = 1; i < N; i++)
	{
		pivot = max_in_column(U, i, i);
		P = make_pivot(P, i, pivot);
		M = make_pivot(M, 0, 0);
		U = matrix_multiplication(P, U);
		
		// Gaussian Elimination
		for (int j = i+1; j < N; j++)
		{
			double val1 = U.get_value(j, i);
			double val2 = U.get_value(i, i);
			scaling_factor = val1 / val2;
			M.set_value(j, i, -scaling_factor);
		}

		U = matrix_multiplication(M, U);
		M = matrix_multiplication(M, P);
		L = matrix_multiplication(M, L);
		A = matrix_multiplication(P, A);		
	}
	// Set P = A. This is the total pivot matrix
	for (int i = 0; i < n * n; i++)
	{
		P.body[i] = A.body[i];
	}

	// Compute the inverse of the total pivot matrix to organize the inverse Lower triangular matrix L. set to A
	matrix B = invert_pivot(A);
	A = matrix_multiplication(L, B);
	B.del_self();

	// Compute the inverse of A. Set to L
	M = make_pivot(M, 0, 0);
	for (int j = 1; j < N; j++)
	{
		for (int i = 1; i < N; i++)
		{
			c.set_value(i, 1, M.get_value(i, j));
		}
		matrix y = forward_substitution(A, c);

		//Fill in A inverse (matrix L)
		for (int i = 1; i < N; i++)
		{
			L.set_value(i, j, y.get_value(i, 1));
		}
		y.del_self();
	}

	matrix ans(n, n * 3);


	//Combine L, U, and P matrices (in that column order) into one matrix to return. Would like to return them separately, don't know how yet.
	for (int i = 1; i < N; i++)
	{
		for (int j = 1; j < N; j++)
		{
			ans.set_value(i, j, L.get_value(i, j));
			ans.set_value(i, j + n, U.get_value(i, j));
			ans.set_value(i, j + 2*n, P.get_value(i, j));
		}
	}

	//Cleanup Extra Matrices
	A.del_self();
	M.del_self();
	c.del_self();
	L.del_self();
	U.del_self();
	P.del_self();

	return ans;
}

void print_matrix(matrix& mat) {
	for (int i = 1; i < mat.nrows+1; i++)
	{
		for (int j = 1; j < mat.ncols + 1; j++)
		{
			double val = mat.get_value(i, j);
			std::cout << val << "    ";
		}
		std::cout << "\n" << "-------------------------------------" << "\n";
	}

	std::cout << "\n" << "\n" << "\n";
}