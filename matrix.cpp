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

matrix LU(matrix& mat) {
	//mat should be a square matrix (size MxM)
	int pivot;
	int n = mat.nrows;
	int N = n + 1;
	double scaling_factor = 1;
	double val = 0;

	matrix P(n, 1);
	matrix L(n, n);
	matrix U(n, n);

	for (size_t i = 0; i < n*n; i++)
	{
		U.body[i] = mat.body[i];
	}

	for (size_t i = 1; i < N; i++)
	{
		pivot = max_in_column(U, i, i);
		P.set_value(pivot, 1, i);

		//Pivot if Neccessary
		if (pivot != i)
		{
			for (size_t j = 1; j < N; j++)
			{
				val = U.get_value(i, j);
				U.set_value(i, j, U.get_value(pivot, j));
				U.set_value(pivot, j, val);
			}
		}
		
		L.set_value(i, i, 1);
		if (i>0)
		{	// Gaussian Elimination
			for (size_t j = i+1; j < N; j++)
			{
				double val1 = U.get_value(j, i);
				double val2 = U.get_value(i, i);
				scaling_factor = val1 / val2;
				L.set_value(j, i, scaling_factor);
				for (size_t k = 1; k < N; k++)
				{
					double val1 = U.get_value(j, k);
					double val2 = U.get_value(i, k);
					val = val1 - scaling_factor * val2;
					U.set_value(j, k, val);
				}
			}
		}
	}
	return matrix_multiplication(L, U);
}