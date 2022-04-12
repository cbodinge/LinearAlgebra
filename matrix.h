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
	double get_value(int row, int col) const;
	void print_matrix() const;
	void fill(double value);
	void del_self();

	void swap_rows(int row1, int row2);
	void identity();

	matrix transpose();
	matrix LU();
	matrix invert_pivot();

	matrix(int nrows=1, int ncols=1) 
	{
		matrix::nrows = nrows;
		matrix::ncols = ncols;
		int size = nrows * ncols;
		body = new double[size];

		fill(0);		
	}
};

matrix addition(matrix& mat1, matrix& mat2);

matrix matrix_multiplication(matrix& mat1, matrix& mat2);

matrix scalar_multiplication(matrix& mat, double scalar);

int max_in_column(matrix& mat, int column, int start);

matrix LU(matrix& mat);

matrix forward_substitution(matrix& lower_trianglular, matrix& b);
