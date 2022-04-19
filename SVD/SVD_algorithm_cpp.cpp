
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string.h>
#include <chrono>  
#include "jacobi.h"

using namespace std;
using namespace std::chrono;


int rows;
int cols;

/* Reads a file line by line - number by number into a vector, saves input dimension (#rows and #colomns) and returns a vector*/
vector<double> getVectorFromFile(string filename)
{
	vector<double> result;
	ifstream input(filename);
	string s;
	rows = 0;
	while (getline(input, s))//read line
	{
		cols = 0;
		istringstream iss(s);
		string num;

		while (getline(iss, num, ','))// read numbers
		{
			result.push_back(stod(num));
			cols++;
		}
		rows++;
	}
	return result;
}

/* Converts vector to 2D-array (matrix), in dimension nxm, and returns it*/
 double** getMatrix(vector<double>& vector, int n, int m)
{
	// initialize 2D-array
	double** M = new double* [n];
	for (int i = 0; i < n; ++i)
		M[i] = new double[m];
	int k = 0;
	// copy vector elemnts to array
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			M[i][j] = vector[k++];
	return M;
}
/* Gets matrix (nxm), returns transpose of matrix (mxn) */
double** transpose(double** M, int n, int m)
{
	// initialize
	double** Mtranspose = new double* [m];
	for (int i = 0; i < m; ++i)
		Mtranspose[i] = new double[n];
	// calc transpose(matrix)
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j) 
			Mtranspose[j][i] = M[i][j];
	
	return Mtranspose;
}
/* Gets 2 matrices: M1(n1xm1) and M2(n2xm2), calculates dot product of M1 and M2 and returns M (n1xm2)  */
double** matricesMultiplication(double** M1, int n1, int m1, double** M2, int n2, int m2)
{
	// rule: n1xm1 mul n2xm2 only if m1==m2
	if (m1 != n2)
	{
		std::cerr << "Can't multiply beetween matrices.\n";
		exit(1);
	}
	// rule: n1xm1 mul n2xm2 = n1xm2
	double** Mresult = new double* [n1];
	for (int i = 0; i < n1; ++i)
		Mresult[i] = new double[m2] {0.0};
	// calc M1 dot M2
	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < m2; j++) {
			for (int k = 0; k < m1; k++)
				Mresult[i][j] += M1[i][k] * M2[k][j];
		}
	}
	return Mresult;
}

/* Gets matrix(nxm) and prints it row by row */
 void printMatrix(double** M, int n, int m, string matrixName)
{
	cout << ">----------------- Print Matrix " << matrixName << " -----------------<" << endl << endl;
	for (int i = 0; i < n; ++i) {
		cout << "[";
		for (int j = 0; j < m; ++j)
			cout << "	" << M[i][j];
		cout << "	]" << endl;
	}
	cout << endl << endl;
}
/* Gets matrix(nxm) and writes it to a file*/
 void writeMatrixToFile(double** M, int n, int m, string matrixName, ofstream& file)
{
	file << ">----------------- Matrix " << matrixName << " -----------------<" << endl << endl;
	for (int i = 0; i < n; ++i) {
		file << "[";
		for (int j = 0; j < m; ++j)
			file << "	" << M[i][j];
		file << "	]" << endl;
	}
	file << endl << endl;
}

/* Gets jacobi object and returns diagonal matrix S with eigenvalues in square root */
double** getS(Jacobi* obj)
{
	int n = obj->N;
	double** S = new double* [n];
	// fill S woth zeros
	for (int i = 0; i < n; i++)
		S[i] = new double[n] {0.0};
	int k = 0;
	for (int i = 0; i < n; i++)
		if (obj->eigenvalues[k] >= 0) // avoid from negative eigen values (zero instead)
			S[i][i] = sqrt(obj->eigenvalues[k++]);
	return S;
}
/*	Calculates pseudo inverse : divide 1 by each non zero value and puts zero for zero  value */
double** getPseudoInverse(double** M, int n)
{
	double** PsInvM = new double* [n];
	for (int i = 0; i < n; i++)
		PsInvM[i] = new double[n] {0.0};
	for (int i = 0; i < n; i++) 
		if (M[i][i] != 0) // avoid from deviding by zero
			PsInvM[i][i] = 1 / M[i][i];
	
	return PsInvM;
}
/* Swaps values between n1 and n2 */
inline void swap(double& n1, double& n2)
{
	double temp = n1;
	n1 = n2;
	n2 = temp;
}

// Gets an array in size n and sorts it with a "Bubble sort" in DESCENDING order
int* bubbleSort(double* arr, int n)
{
	int i, j;
	int* indexes = new int[n];
	for (i = 0; i < n; i++)
		indexes[i] = i; // indexes before sort

	for (i = 0; i < n - 1; i++)
		// Last i elements are already in place
		for (j = 0; j < n - i - 1; j++)
			if (arr[j] < arr[j + 1]) {
				swap(arr[j], arr[j + 1]);
				swap(indexes[j], indexes[j + 1]); // save indexes for eigen vectors sorting (explanation in line 242)
			}
	return indexes;
}
/* Gets Jacobi obj and order array, and returns matrix V with eigenvectors in the order of eigenvalues in S (explanation in line 242) */
double** getV(Jacobi* obj, int* order) {
	int n = obj->N;
	// initialize
	double** V = new double* [n];
	for (int i = 0; i < n; ++i)
		V[i] = new double[n];
	for (int j = 0; j < n; j++) // run on columns
		for (int i = 0; i < n; i++) // run on rows 
			V[i][j] = obj->eigenvectors[i][order[j]]; // V[i][j]=eigenvectors[i][suitable column]
	return V;

}
// if for data_8 equals()=false - we will use thes method for round
inline double round_d(double n, int decimals)
{
	int x = pow(10, decimals); 
	double value = (int)(n * x + .5);
	return (double)value / x;
}
/* Gets 2 matrices: M1(n1xm1) and M2(n2xm2), returns true if equals, else returns false */
inline bool equals(double** M1, int n1, int m1, double** M2, int n2, int m2)
{
	int i, j;
	if (n1 != n2 || m1 != m2)
		return false;
	for (i = 0; i < n1; i++)
		for (j = 0; j < m1; j++)
			if (round_d(M1[i][j],5) != round_d(M2[i][j],5))
				return false;
	return true;
}

int	 main()
{
	//save start time
	auto start = high_resolution_clock::now();

	string file_name = "data_8";

	// read file and insert comma sepereted elemetnts into a vector, update rows and cols (global vars)
	vector<double> input = getVectorFromFile(file_name+".txt");

	// get n x m matrix from n + m vector 
	int a_rows = rows, a_cols = cols;
	double** A = getMatrix(input, a_rows, a_cols);
	//printMatrix(A, rows, cols, "A");


	// -------------------------------- START CALC SVD(M) -------------------------------- //

	cout << "Calculate SVD(A)..." << endl;

	

	// calculate transpose of M
	int at_rows = a_cols, at_cols = a_rows;
	double** At = transpose(A, a_rows, a_cols);
	//printMatrix(At, at_rows, at_cols, "T(At)");

	// Calculate Transpose(M)*M
	double** AtA = matricesMultiplication(At, at_rows, at_cols, A, a_rows, a_cols);
	int ata_size = a_cols;
	//printMatrix(AtA, ata_size, ata_size, "T(A)*A");
	
	//calculate eigen vectors and eigen values by using Jacobi alg on T(M)*M
	Jacobi* jacobiAtA = new Jacobi(AtA, ata_size, ata_size);
	jacobiAtA->calculateEigensByJacobiAlgo();
	//printMatrix(AtA, ata_size, ata_size, "T(A)*A");

	int num_of_eigens = jacobiAtA->N;

	//cout << "eigen values before sort:" << endl;
	//for (int i = 0; i < num_of_eigens; i++)
	//	cout << jacobiAtA->eigenValues[i] << ",	";
	//cout << endl << endl;


	int* order = bubbleSort(jacobiAtA->eigenvalues, num_of_eigens);
	// --->>> Explanation: order includes in each cell what the number of column (eigen vector) that need to be in this place
	//		  for example:	if oreder look like:	cell:	0	1	2	3	so we know that V[:][0] need to be eigenVector[:][3]
	//												val:	3	1	0	2


	//cout << "eigen values after sort:" << endl;
	//for (int i = 0; i < num_of_eigens; i++)
	//	cout << jacobiAtA->eigenValues[i] << ",	";
	//cout << endl << endl;

	//for (int i = 0; i < num_of_eigens; i++)
	//	cout << order[i] << ",	";
	//cout << endl << endl;

	// S
	double** S = getS(jacobiAtA);
	int s_size = num_of_eigens;
	//printMatrix(S, s_size, s_size, "S");

	// V
	int v_size = num_of_eigens;
	double** V = getV(jacobiAtA, order);
	//printMatrix(V, v_size, v_size, "V");

	// pseudo inverse (S) 
	double** s_pseudoinverted = getPseudoInverse(S, s_size);
	//printMatrix(s_pseudoinverted, s_size, s_size, "s_pseudoinverted");

	// U = A(nxm) * V(mxm) * pseudo_inv(S)(mxm)
	double** U = matricesMultiplication(matricesMultiplication(A, a_rows, a_cols, V, v_size, v_size), a_rows, v_size, s_pseudoinverted, s_size, s_size);
	int u_rows = a_rows, u_cols = s_size;
	//printMatrix(U, u_rows, u_cols, "U");

	

	// -------------------------------- FINISH CALC SVD(M) -------------------------------- //

	//printMatrix(S, s_size, s_size, "S");
	//printMatrix(V, v_size, v_size, "V");
	//printMatrix(U, u_rows, u_cols, "U");

	// chack result

	// Vt = T(V)
	double** Vt = transpose(V, v_size, v_size);
	//printMatrix(Vt, v_size, v_size, "Vt");

	// USVt = U * S * T(V)
	double** USVt = matricesMultiplication(matricesMultiplication(U, u_rows, u_cols, S, s_size, s_size), u_rows, s_size, Vt, v_size, v_size);
	//printMatrix(USVt, u_rows, v_size, "USVt");

	string result = "";

	//	 check if U*S*T(V)~= M
	if (equals(A, a_rows, a_cols, USVt, u_rows, v_size))
		result = "U*S*T(V) ~= A";
	else
		result = "U*S*T(V) != A";
	cout << result << endl;
	
	

	// write results to output file
	ofstream output_file;
	output_file.open(file_name.append("_output.txt"));
	writeMatrixToFile(A, a_rows, a_cols, "A", output_file);
	writeMatrixToFile(V, v_size, v_size, "V", output_file);
	writeMatrixToFile(S, s_size, s_size, "S", output_file);
	writeMatrixToFile(U, u_rows, u_cols, "U", output_file);
	writeMatrixToFile(USVt, u_rows, v_size, "USVt", output_file);
	
	// save a finish time
	auto stop = high_resolution_clock::now();
	
	// calc duration of svd calc
	auto duration_sec = duration_cast<seconds>(stop - start);
	auto duration_millis = duration_cast<milliseconds>(stop - start);
	auto duration_nanos = duration_cast<nanoseconds>(stop - start);
	
	cout << "Duration: " << duration_sec.count() << " sec" << endl << endl;
	cout << "Duration: " << duration_millis.count() << " millis" << endl << endl;
	cout << "Duration: " << duration_nanos.count() << " nanos" << endl << endl;
	
	output_file << "Result: " << result << endl;
	output_file << "Duration: " << duration_sec.count() << " sec" << endl << endl;
	output_file << "Duration: " << duration_millis.count() << " millis" << endl << endl;
	output_file << "Duration: " << duration_nanos.count() << " nano" << endl << endl;
	output_file.close();


	

	jacobiAtA->~Jacobi();
	for (int i = 0; i < s_size; i++)
		delete[] S[i];
	delete[] S;
	for (int i = 0; i < v_size; i++)
		delete[] V[i];
	delete[] V;
	for (int i = 0; i < u_rows; i++)
		delete[] U[i];
	delete[] U;
	for (int i = 0; i < a_rows; i++)
		delete[] A[i];
	delete[] A;
	for (int i = 0; i < at_rows; i++)
		delete[] At[i];
	delete[] At;

}