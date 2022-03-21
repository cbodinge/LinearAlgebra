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

class identity : public matrix {
public:
	identity(int nrows, int ncols) {
		identity::nrows = nrows;
		identity::ncols = ncols;
		int size = nrows * ncols;
		body = new double[size];

		for (size_t i = 0; i < size; i++)
		{
			int row = floor(i / identity::ncols) + 1;
			int col = i % ncols + 1;
			identity::set_value(row, col, 0);
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

void print_matrix(matrix& mat);

matrix make_pivot(matrix& mat, int row1, int row2);

matrix invert_pivot(matrix& mat);

matrix forward_substitution(matrix& lower_trianglular, matrix& b);
