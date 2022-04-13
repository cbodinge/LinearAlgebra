#pragma once
#include <iostream>
#include <vector>

class Matrix
{
public:
	double* body;
	int nrows;
	int ncols;

	void set_value(int row, int col, double value);
	double get_value(int row, int col) const;
	void print_matrix() const;
	void fill(double value);

	void swap_rows(int row1, int row2);
	void identity();

	Matrix transpose();
	Matrix LU();
	Matrix invert_pivot();

	Matrix(int nrows=1, int ncols=1) 
	{
		Matrix::nrows = nrows;
		Matrix::ncols = ncols;
		int size = nrows * ncols;
		body = new double[size];

		fill(0);		
	}

	//Destructor
	~Matrix()
	{ delete[] body; }
};

Matrix addition(Matrix& mat1, Matrix& mat2);

Matrix matrix_multiplication(Matrix& mat1, Matrix& mat2);

Matrix scalar_multiplication(Matrix& mat, double scalar);

int max_in_column(Matrix& mat, int column, int start);

Matrix LU(Matrix& mat);

Matrix forward_substitution(Matrix& lower_trianglular, Matrix& b);
