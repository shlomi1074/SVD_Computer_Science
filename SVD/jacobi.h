/*
The C++ class of the classical Jacobi's method implementation for finding eigenvalues and
eigenvectors of a symetric matrix which represented as 1-dimension array.
*/

#include <string>
#ifndef JACOBI_H
#define JACOBI_H

class Jacobi {

private:
	const double tol;

	double** M = NULL;

	int k;
	int l;

	double findMaxElem();
	void rotate();

	void printMatrix(double** M, int n, int m, std::string matrixName);

public:
	const int N;

	double** eigenvectors = NULL;
	double* eigenvalues = NULL;

	Jacobi(double** input_matrix, int input_size, const double tol = 1.0e-9);
	~Jacobi();

	void calculateEigensByJacobiAlgo();

};

#endif#pragma once
