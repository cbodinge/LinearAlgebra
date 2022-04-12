#include "matrix.h"


//-----------------------------------------------------------------------------------------------------------
//	Class Methods for matrix
//-----------------------------------------------------------------------------------------------------------

void matrix::set_value(int row, int col, double value)
{
	int position = (row - 1) * (matrix::ncols) + col - 1;
	matrix::body[position] = value;
}

double matrix::get_value(int row, int col) const
{
	int position = (row - 1) * (matrix::ncols)+col - 1;
	double ans = body[position];
	return ans;
}

void matrix::print_matrix() const
{
	for (int i = 1; i < nrows + 1; i++)
	{
		for (int j = 1; j < ncols + 1; j++)
		{
			double val = get_value(i, j);
			val = round(val * 1000) / 1000;
			std::cout << val << "\t";
		}
		std::cout << "\n" << "-----------------------------------------------------------------" << "\n";
	}

	std::cout << "\n" << "\n" << "\n";
}


matrix matrix::transpose()
{
	matrix ans(ncols, nrows);
	double val;
	for (int i = 1; i < nrows + 1; i++)
	{
		for (int j = 1; j < ncols + 1; j++)
		{
			val = get_value(i, j);
			ans.set_value(j, i, val);
		}
	}
	return ans;
}

void matrix::identity()
{
	for (int i = 1; i < nrows + 1; i++)
	{
		for (int j = 1; j < ncols + 1; j++)
		{
			if (i == j)
			{
				set_value(i, j, 1);
			}
			else
			{
				set_value(i, j, 0);
			}
		}
	}
}

void matrix::swap_rows(int row1, int row2)
{
	double val1;
	double val2;
	for (int j = 1; j < ncols + 1; j++)
	{
		val1 = get_value(row1, j);
		val2 = get_value(row2, j);
		set_value(row2, j, val1);
		set_value(row1, j, val2);
	}
}

matrix matrix::invert_pivot()
	// Must be a pivot matrix to return a correct result.
{
	matrix temp(nrows, ncols);
	for (int i = 1; i < nrows + 1; i++)
	{
		for (int j = 1; j < ncols + 1; j++)
		{
			temp.set_value(nrows + 1 - i, ncols + 1 - j, get_value(i, j));
		}
	}

	return temp;
}

matrix matrix::LU() {
	//mat should be a square matrix (size MxM)
	int pivot;
	int n = nrows;
	int N = n + 1;
	int size = n * n;
	double scaling_factor = 1;
	double val = 0;

	matrix complete_pivot_matrix(n, n);
	matrix upper_matrix(n, n);
	matrix lower_matrix(n, n);

	matrix pivot_matrix(n, n);
	matrix row_operations(n, n);

	for (int i = 0; i < n * n; i++)
	{
		upper_matrix.body[i] = body[i];
	}

	lower_matrix.identity();
	complete_pivot_matrix.identity();
	for (int i = 1; i < N; i++)
	{
		pivot = max_in_column(upper_matrix, i, i);
		pivot_matrix.identity();
		pivot_matrix.swap_rows(i, pivot);
		upper_matrix.swap_rows(i, pivot);
		lower_matrix.swap_rows(i, pivot);
		complete_pivot_matrix.swap_rows(i, pivot);

		// Gaussian Elimination
		row_operations.identity();
		for (int j = i + 1; j < N; j++)
		{
			double val1 = upper_matrix.get_value(j, i);
			double val2 = upper_matrix.get_value(i, i);
			scaling_factor = val1 / val2;
			row_operations.set_value(j, i, -scaling_factor);
		}

		// Update Upper Matrix with Row Operations
		matrix temp = matrix_multiplication(row_operations, upper_matrix);
		for (int i = 0; i < size; i++)
		{
			upper_matrix.body[i] = temp.body[i];
		}

		// Update Upper Matrix with Row Operations
		temp = matrix_multiplication(row_operations, lower_matrix);
		for (int i = 0; i < size; i++)
		{
			lower_matrix.body[i] = temp.body[i];
		}

	}

	// Compute the inverse of the total pivot matrix to organize the inverse Lower triangular matrix L.
	matrix inverted_pivot = complete_pivot_matrix.invert_pivot();
	matrix inverted_lower_matrix = matrix_multiplication(lower_matrix, inverted_pivot);

	matrix temp(n, 1);
	matrix iden(n, n);
	iden.identity();

	for (int j = 1; j < N; j++)
	{
		for (int i = 1; i < N; i++)
		{
			temp.set_value(i, 1, iden.get_value(i, j));
		}
		matrix y = forward_substitution(inverted_lower_matrix, temp);

		// Fill in A inverse (matrix L)
		for (int i = 1; i < N; i++)
		{
			lower_matrix.set_value(i, j, y.get_value(i, 1));
		}
	}
	matrix ans(n, n * 3);


	//Combine L, U, and P matrices (in that column order) into one matrix to return. Would like to return them separately, don't know how yet.
	for (int i = 1; i < N; i++)
	{
		for (int j = 1; j < N; j++)
		{
			ans.set_value(i, j, lower_matrix.get_value(i, j));
			ans.set_value(i, j + n, upper_matrix.get_value(i, j));
			ans.set_value(i, j + 2 * n, complete_pivot_matrix.get_value(i, j));
		}
	}

	return ans;
}

void matrix::fill(double value)
{
	for (int i = 1; i < nrows + 1; i++)
	{
		for (int j = 1; j < ncols + 1; j++)
		{
			set_value(i, j, value);
		}
	}
}


// ----------------------------------------------------------------------------------------------------------
//	Functions Involving the matrix class
// ----------------------------------------------------------------------------------------------------------

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

matrix forward_substitution(matrix& lower_trianglular, matrix& b)
{
	matrix y(b.nrows, 1);
	double inner_sum;
	y.set_value(1, 1, b.get_value(1, 1));

	for (int i = 2; i < y.nrows + 1; i++)
	{
		inner_sum = 0;
		for (int j = 0; j < i; j++)
		{
			double val1 = lower_trianglular.get_value(i, j);
			double val2 = y.get_value(j, 1);
			inner_sum = inner_sum + val1 * val2;
			int a = 1;
		}
		y.set_value(i, 1, b.get_value(i, 1)-inner_sum);
	}

	return y;
}
