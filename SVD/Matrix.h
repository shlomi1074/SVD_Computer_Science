#pragma once

namespace LA
{
	class Matrix
	{
	public:
		Matrix(int nRows, int nCols, bool isRandom);
		Matrix(const Matrix& inputMatrix);
		Matrix Transpose();
		Matrix PseudoInverse();
		 
		void CalculateEigensJacobi();
		void SetElement(int row, int col, float elementValue);
		void SetElement(int index, float elementValue);

		void PrintMatrix();
		bool Compare(const Matrix& matrix1, float tolerance);
		bool IsSquare();
		float GetElement(int row, int col);
		//int Rank();
		//~Matrix();

		bool operator== (const Matrix& rhs);
		Matrix operator= (const Matrix& rhs);
		friend Matrix operator+ (const Matrix& lhs, const Matrix& rhs);
		friend Matrix operator- (const Matrix& lhs, const Matrix& rhs);
		friend Matrix operator* (const Matrix& lhs, const Matrix& rhs);


	private:
		int calculateIndex(int row, int col);
		float m_maxValue;
		int m_maxValueRowIndex;
		int m_maxValueColIndex;
		float* m_matrix;
		Matrix* m_eigenVectors = nullptr;
		Matrix* m_eigenValues = nullptr;
		int m_nRows, m_nCols, m_nElements;

	};
}
