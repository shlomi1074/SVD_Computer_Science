#pragma once
#include <vector>
#include <iostream>

using namespace std;

namespace LA
{
	class Matrix
	{
	public:
		Matrix(int nRows, int nCols, bool isRandom);
		Matrix(int nRows, int nCols, vector<long double>* inputData);
		Matrix(const Matrix& inputMatrix);

		Matrix Transpose();
		Matrix PseudoInverse();
		void SVD(long double tolerance);
		inline Matrix GetVt() { return *m_eigenVectorsTranspose; }
		inline Matrix GetV() { return *m_eigenVectors; }
		inline Matrix GetS() { return *m_eigenValues; }
		inline Matrix GetSi() { return *m_eigenValuesInverse; }
		inline Matrix GetU() { return (*this) * (*m_eigenVectors) * (*m_eigenValuesInverse); }

		void CalculateEigensJacobi(long double tolerance);
		void SetElement(int row, int col, long double elementValue);
		void SetElement(int index, long double elementValue);
		void WriteMatrixToFile(string fileName, string matrixName);
		void PrintMatrix();
		bool Compare(const Matrix& matrix1, long double tolerance);
		bool IsSquare();
		long double GetElement(int row, int col) const;
		long double findMaxElem();

		inline int GetColNum() { return m_nCols; }
		inline int GetRowNum() { return m_nRows; }

		bool operator== (const Matrix& rhs);
		Matrix operator= (const Matrix& rhs);
		friend Matrix operator+ (const Matrix& lhs, const Matrix& rhs);
		friend Matrix operator- (const Matrix& lhs, const Matrix& rhs);
		friend Matrix operator* (const Matrix& lhs, const Matrix& rhs);

	private:
		int calculateIndex(int row, int col) const;
		void jacobiEigensRotate();
		void arrangeEigenVectors(vector<pair<long double, int>> eigenValues);
		long double m_maxValue;
		int m_maxValueRowIndex;
		int m_maxValueColIndex;
		long double* m_matrix;
		Matrix* m_eigenVectors = nullptr;
		Matrix* m_eigenVectorsTranspose = nullptr;
		Matrix* m_eigenValues = nullptr;
		Matrix* m_eigenValuesInverse = nullptr;

		int m_nRows, m_nCols, m_nElements;

	};
}
