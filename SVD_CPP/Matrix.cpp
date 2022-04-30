#include "Matrix.h"
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <math.h>

LA::Matrix::Matrix(int nRows, int nCols, bool isRandom):
	m_nRows(nRows),
	m_nCols(nCols),
	m_nElements(m_nRows* m_nCols)
{
	m_maxValueRowIndex = -1;
	m_maxValueColIndex = -1;
	m_matrix = new double[m_nElements];
	if (isRandom) {
		srand(time(NULL));
		for (int i = 0; i < m_nElements; i++)
		{
			SetElement(i, rand() % 1000 + 700);
		}
	}
	else {
		for (int i = 0; i < m_nElements; i++)
		{
			SetElement(i, 0);
		}
	}
}

LA::Matrix::Matrix(int nRows, int nCols, vector<double>* inputData):
m_nRows(nRows),
m_nCols(nCols),
m_nElements(m_nRows* m_nCols)
{
	m_matrix = new double[m_nElements];
	for (int i = 0; i < m_nElements; i++)
		m_matrix[i] = (*inputData)[i];

	m_maxValueRowIndex = -1;
	m_maxValueColIndex = -1;
	m_maxValue = std::numeric_limits<double>::min();
	for (int i = 0; i < m_nRows; i++)
	{
		for (int j = 0; j < m_nCols; j++)
		{
			if (i != j && fabs(GetElement(i, j)) > m_maxValue) 
			{
				m_maxValue = fabs(GetElement(i, j));
				m_maxValueRowIndex = i;
				m_maxValueColIndex = j;
			}
		}
	}
}

LA::Matrix::Matrix(const Matrix& inputMatrix) :
	m_nRows(inputMatrix.m_nRows),
	m_nCols(inputMatrix.m_nCols),
	m_nElements(inputMatrix.m_nElements),
	m_maxValue(inputMatrix.m_maxValue)
{
	m_maxValueRowIndex = inputMatrix.m_maxValueRowIndex;
	m_maxValueColIndex = inputMatrix.m_maxValueColIndex;
	m_eigenValues = inputMatrix.m_eigenValues;
	m_eigenVectors = inputMatrix.m_eigenVectors;
	m_eigenValuesInverse = inputMatrix.m_eigenValuesInverse;
	m_eigenVectorsTranspose = inputMatrix.m_eigenVectorsTranspose;

	m_matrix = new double[m_nElements];
	for (int i = 0; i < m_nElements; i++)
		m_matrix[i] = inputMatrix.m_matrix[i];
}

double LA::Matrix::GetElement(int row, int col) const
{
	auto index = LA::Matrix::calculateIndex(row, col);
	return m_matrix[index];
}

void LA::Matrix::CalculateEigensJacobi(double tolerance = 1.0e-8)
{
	//TODO: check if ok
	m_eigenValues = new LA::Matrix(m_nCols, m_nCols, false);
	m_eigenVectors = new LA::Matrix(m_nCols, m_nCols, false);
	m_eigenValuesInverse = new LA::Matrix(m_nCols, m_nCols, false);
	m_eigenVectorsTranspose = new LA::Matrix(m_nCols, m_nCols, false);
	for (int i = 0; i < m_nCols; i++)
	{
		m_eigenVectors->SetElement(i, i, 1);
	}

	int numOfIterations = 2 * m_nCols;
	//const int tol = tolerance;
	for (int i = 0; i < numOfIterations; i++) {
		double aMax = findMaxElem();
		if (m_maxValue < tolerance)
			break;
		jacobiEigensRotate();
	}
	for (int i = 0; i < m_nCols; i++)
		m_eigenValues->SetElement(i, i, GetElement(i, i));
}

void LA::Matrix::jacobiEigensRotate()
{
	//calculations
	double maxRow =  GetElement(m_maxValueRowIndex, m_maxValueRowIndex);
	double maxCol = GetElement(m_maxValueColIndex, m_maxValueColIndex);
	double maxRowCol = GetElement(m_maxValueRowIndex, m_maxValueColIndex);

	double aDiff = maxCol - maxRow;
	double t, phi;
	if (fabs(maxRowCol) < fabs(aDiff) * 1.0e-36)
		t = maxRowCol / aDiff;
	else {
		phi = aDiff / (2.0 * maxRowCol);
		t = 1.0 / (fabs(phi) + sqrt(pow(phi, 2) + 1.0));
		if (phi < 0.0)
			t = -t;
	}
	double c = 1.0 / sqrt(pow(t, 2) + 1.0);
	double s = t * c;
	double tau = s / (1.0 + c);
	double temp = maxRowCol;
	double tempValue;
	double elem;
	double tempCol;
	double tempRow;
	//updating M (for eigen values)
 	SetElement(m_maxValueRowIndex, m_maxValueColIndex, 0);
	SetElement(m_maxValueRowIndex, m_maxValueRowIndex, maxRow - t * temp);
	SetElement(m_maxValueColIndex, m_maxValueColIndex, maxCol + t * temp);

	for (int i = 0; i < m_maxValueRowIndex; i++) {// Case of i < k
		temp = GetElement(i, m_maxValueRowIndex);
		elem = GetElement(i, m_maxValueColIndex);
		tempValue = temp - s * (elem + tau * temp);
		SetElement(i, m_maxValueRowIndex, tempValue);
		elem = GetElement(i, m_maxValueColIndex);
		tempValue = elem + s * (temp - tau * elem);
		SetElement(i, m_maxValueColIndex, tempValue);
	}
	for (int i = m_maxValueRowIndex + 1; i < m_maxValueColIndex; i++) {// Case of k < i < l
		temp = GetElement(m_maxValueRowIndex, i);
		elem = GetElement(i, m_maxValueColIndex);
		tempValue = temp - s * (elem + tau * temp);
		SetElement(m_maxValueRowIndex, i, tempValue);
		elem = GetElement(i, m_maxValueColIndex);
		tempValue = elem + s * (temp - tau * elem);
		SetElement(i, m_maxValueColIndex, tempValue);
	}
	for (int i = m_maxValueColIndex + 1; i < m_nCols; i++) { // Case of i > l
		temp = GetElement(m_maxValueRowIndex, i);
		elem = GetElement(m_maxValueColIndex, i);
		tempValue = temp - s * (elem + tau * temp);
		SetElement(m_maxValueRowIndex, i, tempValue);
		elem = GetElement(m_maxValueColIndex, i);
		tempValue = elem + s * (temp - tau * elem);
		SetElement(m_maxValueColIndex, i, tempValue);
	}
	//updating eigenvectors
	for (int i = 0; i < m_nCols; i++) {//Update transformation matrix
		tempRow = m_eigenVectors->GetElement(i, m_maxValueRowIndex); 
		tempCol = m_eigenVectors->GetElement(i, m_maxValueColIndex); 

		tempValue = tempRow - s * (tempCol + tau * tempRow);
		m_eigenVectors->SetElement(i, m_maxValueRowIndex, tempValue);

		tempCol = m_eigenVectors->GetElement(i, m_maxValueColIndex);
		tempValue = tempCol + s * (tempRow - tau * tempCol);
		m_eigenVectors->SetElement(i, m_maxValueColIndex, tempValue);
	}
}

void LA::Matrix::arrangeEigenVectors(vector<pair<double, int>> eigenValues)
{
	// initialize
	Matrix result(m_nCols, m_nCols, false);
	for (int j = 0; j < m_nCols; j++) // run on columns
		for (int i = 0; i < m_nRows; i++) { // run on rows
			double elem = m_eigenVectors->GetElement(i, eigenValues[j].second);
			result.SetElement(i, j, elem);
		}
	*m_eigenVectors = result;	
}

void LA::Matrix::SetElement(int row, int col, double elementValue)
{
	auto index = LA::Matrix::calculateIndex(row, col);
	m_matrix[index] = elementValue;
}

void LA::Matrix::SetElement(int index, double elementValue)
{
	m_matrix[index] = elementValue;
}

void LA::Matrix::WriteMatrixToFile(string fileName, string matrixName)
{
	std::ofstream output_file;
	output_file.open(fileName.append("_output.txt"), std::ios_base::app);
	output_file << ">----------------- Matrix " << matrixName << " -----------------<" << endl << endl;
	for (int i = 0; i < m_nRows; ++i) {
		output_file << "[";
		for (int j = 0; j < m_nCols; ++j)
			output_file << "	" << GetElement(i, j);
		output_file << "	]" << endl;
	}
	output_file << endl << endl;
	output_file.close();
}

int LA::Matrix::calculateIndex(int row, int col) const
{
	return (row * m_nCols) + col;
}

LA::Matrix LA::Matrix::Transpose()
{
	Matrix TMatrix(m_nCols, m_nRows, false);

	for (int i = 0; i < m_nRows; ++i)
	{
		for (int j = 0; j < m_nCols; ++j)
		{
			TMatrix.SetElement(j, i, this->GetElement(i, j));
		}
	}
	return TMatrix;
}

/// <summary>
/// Works only on diagonal matrices 
/// </summary>
/// <returns></returns>
LA::Matrix LA::Matrix::PseudoInverse()
{
	if (!IsSquare()) {
		throw std::invalid_argument("[PseudoInverse] matrix has to be square");
	}
	Matrix resMatrix(m_nCols, m_nRows, false);

	for (int i = 0; i < m_nCols; i++)
	{
		auto diagValue = GetElement(i, i);
		if (diagValue != 0)
		{
				resMatrix.SetElement(i, i, 1 / diagValue);
		}
	}
	return resMatrix;
}

void LA::Matrix::SVD(double tolerance = 1.0e-8)
{
	m_eigenValues = new LA::Matrix(m_nCols, m_nCols, false);
	m_eigenVectors = new LA::Matrix(m_nCols, m_nCols, false);
	m_eigenValuesInverse = new LA::Matrix(m_nCols, m_nCols, false);
	m_eigenVectorsTranspose = new LA::Matrix(m_nCols, m_nCols, false);

	auto At = Transpose();
	auto AtA = At * (*this);
	AtA.CalculateEigensJacobi(tolerance);
	LA::Matrix S(AtA);

	vector<pair<double, int>> eigenValues;

	for (int i = 0; i < m_nCols; i++)
	{
		eigenValues.push_back(pair<double, int>(S.m_eigenValues->GetElement(i, i), i));
	}
	sort(eigenValues.begin(), eigenValues.end(), greater<pair<double, int>>());
	S.arrangeEigenVectors(eigenValues);
	for (int i = 0; i < m_nCols; i++)
	{
		S.m_eigenValues->SetElement(i, i, sqrt(eigenValues[i].first));
	}
	*m_eigenValuesInverse = S.m_eigenValues->PseudoInverse();
	*m_eigenVectorsTranspose = S.m_eigenVectors->Transpose();
	*m_eigenValues = *S.m_eigenValues;
	*m_eigenVectors = *S.m_eigenVectors;
}

double LA::Matrix::findMaxElem() {
	double aMax = 0.0;
	for (int i = 0; i < m_nCols - 1; i++)
		for (int j = i + 1; j < m_nCols; j++)
			if (aMax <= fabs(GetElement(i, j))) {
				aMax = fabs(GetElement(i, j));
				//save location (i,j) of max element in matrix
				m_maxValueRowIndex = i;
				m_maxValueColIndex = j;
				m_maxValue = aMax;
			}
	return aMax;
}

void LA::Matrix::PrintMatrix()
{
	int nRows = m_nRows;
	int nCols = m_nCols;
	for (int row = 0; row < nRows; ++row)
	{
		for (int col = 0; col < nCols; ++col)
		{
			std::cout << std::fixed << std::setprecision(3) << this->GetElement(row, col) << "  ";
		}
		std::cout << std::endl;
	}
}

bool LA::Matrix::Compare(const Matrix& matrix1, double tolerance)
{
	int numRows1 = matrix1.m_nRows;
	int numCols1 = matrix1.m_nCols;
	if ((numRows1 != m_nRows) || (numCols1 != m_nCols))
		return false;

	double cumulativeSum = 0.0;
	for (int i = 0; i < m_nElements; ++i)
	{
		double element1 = matrix1.m_matrix[i];
		double element2 = m_matrix[i];
		double diff = sqrt(((element1 - element2) * (element1 - element2)));

		if (diff > tolerance)
		{
			return false;
		}
	}
	return true;
}

bool LA::Matrix::IsSquare()
{
	return m_nCols == m_nRows;
}

bool LA::Matrix::operator==(const Matrix& rhs)
{
	if ((this->m_nRows != rhs.m_nRows) && (this->m_nCols != rhs.m_nCols))
		return false;

	bool flag = true;
	for (int i = 0; i < this->m_nElements; ++i)
	{
		auto isCloseEnough = fabs(this->m_matrix[i] - rhs.m_matrix[i]) < 1e-8;
		if (!isCloseEnough)
			flag = false;
	}
	return flag;
}

/// <summary>
/// A = B;
/// Destroying left matrix and create it again
/// </summary>
/// <param name="rhs"></param>
/// <returns></returns>
LA::Matrix LA::Matrix::operator=(const Matrix& rhs)
{
	if (this != &rhs)
	{
		m_nRows = rhs.m_nRows;
		m_nCols = rhs.m_nCols;
		m_nElements = rhs.m_nElements;
		m_maxValue = rhs.m_maxValue;

		if (m_matrix)
			delete[] m_matrix;

		m_matrix = new double[m_nElements];
		for (int i = 0; i < m_nElements; i++)
			m_matrix[i] = rhs.m_matrix[i];

		return *this;
	}

}

LA::Matrix LA::operator*(const Matrix& lhs, const Matrix& rhs)
{
	int r_numRows = rhs.m_nRows;
	int r_numCols = rhs.m_nCols;
	int l_numRows = lhs.m_nRows;
	int l_numCols = lhs.m_nCols;

	if (l_numCols == r_numRows)
	{
		Matrix result(lhs.m_nRows, rhs.m_nCols, false);

			for (int i = 0; i < l_numRows; i++) {
				for (int j = 0; j < r_numCols; j++) {
					for (int k = 0; k < l_numCols; k++) {
						double elementResult = result.GetElement(i, j) + lhs.GetElement(i, k) * rhs.GetElement(k, j);
						result.SetElement(i, j, elementResult);
					}
				}
			}				
		return result;
	}
	else
	{
		Matrix result(1, 1, false);
		result.SetElement(0, 0);
		return result;
	}
}
