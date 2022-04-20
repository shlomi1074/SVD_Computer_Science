#pragma once
#include "Matrix.h"
#include <vector>
#include <algorithm>

namespace LA
{
	class Utilites
	{
	public:
		static vector<pair<float, int>> bubbleSort(vector<pair<float, int>> eigenValues, int n)
		{
			sort(eigenValues.begin(), eigenValues.end());
			return 	eigenValues;
		}
	private:


	};
}