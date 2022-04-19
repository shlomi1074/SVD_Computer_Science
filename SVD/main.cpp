//
//#include <iostream>
//#include "Matrix.h"
//
//int main()
//{ 
//	LA::Matrix mat(2, 2, false);
//	mat.SetElement(0, 0, 6.3245);
//	mat.SetElement(1, 1, 3.1622);
//
//	LA::Matrix mat2(10, 10, false);
//	mat2.SetElement(1, 1, 2);
//	try
//	{
//		LA::Matrix inv3 = mat.PseudoInverse();
//		inv3.PrintMatrix();
//
//	}
//	catch (const std::exception& e)
//	{
//		std::cout << e.what() << std::endl;
//	}
//
//	//auto mat3 = mat + mat2;
//	//mat2.SetElement(1, 2, 2);
//	//auto temp = mat == mat2;
//	//mat.PrintMatrix();
//	//mat.PrintMatrix();
//	//mat3.PrintMatrix();
//	//std::cout << std::endl;
//	//mat2.PrintMatrix();
//
//	//mat.Compare(mat2, 0.001);
//}
