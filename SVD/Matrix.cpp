#include "Matrix.h"
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <iomanip>

LA::Matrix::Matrix(int nRows, int nCols, bool isRandom):
	m_nRows(nRows),
	m_nCols(nCols),
	m_nElements(m_nRows* m_nCols)
{
	m_matrix = new float[m_nElements];

	if (isRandom) {
		srand(time(NULL));
		for (int i = 0; i < m_nElements; i++)
		{
			SetElement(i, rand() % 100);
		}
	}
	else {
		for (int i = 0; i < m_nElements; i++)
		{
			SetElement(i, 0);
		}
	}
}

LA::Matrix::Matrix(const Matrix& inputMatrix) :
	m_nRows(inputMatrix.m_nRows),
	m_nCols(inputMatrix.m_nCols),
	m_nElements(inputMatrix.m_nElements),
	m_maxValue(inputMatrix.m_maxValue)
{


	m_matrix = new float[m_nElements];
	for (int i = 0; i < m_nElements; i++)
		m_matrix[i] = inputMatrix.m_matrix[i];
}

float LA::Matrix::GetElement(int row, int col)
{
	auto index = LA::Matrix::calculateIndex(row, col);
	return m_matrix[index];
}

void LA::Matrix::CalculateEigensJacobi()
{
	//TODO: check if ok
	m_eigenValues = new LA::Matrix(m_nCols, m_nCols, false);
	m_eigenVectors = new LA::Matrix(m_nCols, m_nCols, false);

	
}

void LA::Matrix::SetElement(int row, int col, float elementValue)
{
	auto index = LA::Matrix::calculateIndex(row, col);
	m_matrix[index] = elementValue;

	if (row != col && m_maxValue < m_matrix[index])
	{
		m_maxValue = m_matrix[index];
		m_maxValueRowIndex = row;
		m_maxValueColIndex = col;
	}
}

void LA::Matrix::SetElement(int index, float elementValue)
{
	m_matrix[index] = elementValue;

	int row = index / m_nCols;
	int col = index % m_nCols;
	if (row != col && m_maxValue < m_matrix[index])
	{
		m_maxValue = m_matrix[index];
		m_maxValueRowIndex = row;
		m_maxValueColIndex = col;
	}
}

int LA::Matrix::calculateIndex(int row, int col)
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

bool LA::Matrix::Compare(const Matrix& matrix1, float tolerance)
{
	int numRows1 = matrix1.m_nRows;
	int numCols1 = matrix1.m_nCols;
	if ((numRows1 != m_nRows) || (numCols1 != m_nCols))
		return false;

	double cumulativeSum = 0.0;
	for (int i = 0; i < m_nElements; ++i)
	{
		float element1 = matrix1.m_matrix[i];
		float element2 = m_matrix[i];
		float diff = sqrt(((element1 - element2) * (element1 - element2)));

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

		m_matrix = new float[m_nElements];
		for (int i = 0; i < m_nElements; i++)
			m_matrix[i] = rhs.m_matrix[i];

		return *this;
	}

}

LA::Matrix LA::operator+(const Matrix& lhs, const Matrix& rhs)
{
	int numRows = lhs.m_nRows;
	int numCols = lhs.m_nCols;

	Matrix result(numRows, numCols, false);

	for (int i = 0; i < result.m_nElements; i++)
		result.SetElement(i, lhs.m_matrix[i] + rhs.m_matrix[i]);

	return result;
}

LA::Matrix LA::operator-(const Matrix& lhs, const Matrix& rhs)
{
	int numRows = lhs.m_nRows;
	int numCols = lhs.m_nCols;

	Matrix result(numRows, numCols, false);

	for (int i = 0; i < result.m_nElements; i++)
		result.SetElement(i, lhs.m_matrix[i] - rhs.m_matrix[i]);

	return result;
}

LA::Matrix LA::operator*(const Matrix& lhs, const Matrix& rhs)
{
	int numRows = lhs.m_nRows;
	int numCols = lhs.m_nCols;

	Matrix result(numRows, numCols, false);

	for (int i = 0; i < result.m_nElements; i++)
		result.SetElement(i, lhs.m_matrix[i] * rhs.m_matrix[i]);

	return result;
}
