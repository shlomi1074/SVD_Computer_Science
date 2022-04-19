#include <cmath>
#include <string>
#include <cstdlib>
#include <sstream>
#include <chrono> 
#include <iostream>
#include <math.h>
#include <string.h>
#include "jacobi.h"
using namespace std;


Jacobi::Jacobi(double** input_matrix, int input_size, const double tol) : tol(tol), N(input_size) {
	std::cout << "TOL: " << tol << std::endl;
	M = input_matrix;
	//Initialize eigen vectors 2d array to be Identity matrix
	eigenvectors = new double* [N];
	for (int i = 0; i < N; i++)
		eigenvectors[i] = new double[N] {0.0};
	for (int i = 0; i < N; i++)
		eigenvectors[i][i] = 1.0;

	//initialize eigen values array
	eigenvalues = new double[N];

}

Jacobi::~Jacobi() {
	for (int i = 0; i < N; i++)
		delete[] M[i];
	delete[] M;
	for (int i = 0; i < N; i++)
		delete[] eigenvectors[i];
	delete[] eigenvectors;

	delete[] eigenvalues;
}

void Jacobi::calculateEigensByJacobiAlgo() {

	//find eignes
	int maxRot = 2 * N;
	double aMax;
	for (int i = 0; i < maxRot; i++) {
		aMax = findMaxElem();
		if (aMax < tol)
			break;
		rotate();
	}
	//fills array with eigen values that calculated during rotations
	for (int i = 0; i < N; i++)
		eigenvalues[i] = M[i][i];
	//printMatrix(M, N, N, "M");




}
/* Find max element in matrix (exclude diagonal)*/
double Jacobi::findMaxElem() {
	double aMax = 0.0;
	for (int i = 0; i < N - 1; i++)
		for (int j = i + 1; j < N; j++)
			if (aMax <= fabs(M[i][j])) {
				aMax = fabs(M[i][j]);
				//save location (i,j) of max element in matrix
				k = i;
				l = j;
			}
	return aMax;
}

void Jacobi::rotate() {

	//calculations
	double aDiff = M[l][l] - M[k][k];
	double t, phi;
	if (fabs(M[k][l]) < fabs(aDiff) * 1.0e-36)
		t = M[k][l]/aDiff;
	else {
		phi = aDiff / (2.0 * M[k][l]);
		t = 1.0 / (fabs(phi) + sqrt(pow(phi, 2) + 1.0));
		if (phi < 0.0)
			t = -t;
	}
	double c = 1.0 / sqrt(pow(t, 2) + 1.0);
	double s = t * c;
	double tau = s / (1.0 + c);
	double temp = M[k][l];

	//updating M (for eigen values)
	M[k][l] = 0.0;
	M[k][k] = M[k][k] - t * temp;
	M[l][l] = M[l][l] + t * temp;
	for (int i = 0; i < k; i++) {// Case of i < k
		temp = M[i][k];
		M[i][k] = temp - s * (M[i][l] + tau * temp);
		M[i][l] = M[i][l] + s * (temp - tau * M[i][l]);
	}
	for (int i = k + 1; i < l; i++) {// Case of k < i < l
		temp = M[k][i];
		M[k][i] = temp - s * (M[i][l] + tau * M[k][i]);
		M[i][l] = M[i][l] + s * (temp - tau * M[i][l]);
	}
	for (int i = l + 1; i < N; i++) { // Case of i > l
		temp = M[k][i];
		M[k][i] = temp - s * (M[l][i] + tau * temp);
		M[l][i] = M[l][i] + s * (temp - tau * M[l][i]);
	}
	//updating eigenvectors
	for (int i = 0; i < N; i++) {//Update transformation matrix
		temp = eigenvectors[i][k];
		eigenvectors[i][k] = temp - s * (eigenvectors[i][l] + tau * eigenvectors[i][k]);
		eigenvectors[i][l] = eigenvectors[i][l] + s * (temp - tau * eigenvectors[i][l]);
	}

}
void Jacobi::printMatrix(double** M, int n, int m, string matrixName) {
	cout << ">----------------- Print Matrix " << matrixName << " -----------------<" << endl << endl;
	for (int i = 0; i < n; ++i) {
		cout << "[";
		for (int j = 0; j < m; ++j)
			cout << "	" << M[i][j];
		cout << "	]" << endl;
	}
	cout << endl << endl;
}