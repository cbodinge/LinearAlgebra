#pragma once
#include <iostream>
#include <vector>

class matrix
{
public:
	double* body;
	int nrows;
	int ncols;

	void set_value(int row, int col, double value);
	double get_value(int row, int col);
	void del_self();

	matrix(int nrows, int ncols) {
		matrix::nrows = nrows;
		matrix::ncols = ncols;
		int size = nrows * ncols;
		body = new double[size];

		for (size_t i = 0; i < size; i++)
		{	
			int row = floor(i / matrix::ncols) + 1;
			int col = i % ncols + 1;
			matrix::set_value(row, col, 0);
		}
	}
};

matrix addition(matrix& mat1, matrix& mat2);

matrix matrix_multiplication(matrix& mat1, matrix& mat2);

matrix scalar_multiplication(matrix& mat, double scalar);

matrix transpose(matrix& mat);

matrix solve(matrix& mat1, matrix& mat2);

matrix LU(matrix& mat);

int max_in_column(matrix& mat, int column, int start);

