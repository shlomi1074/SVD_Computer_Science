
#include <iostream>
#include "Matrix.h"
#include <fstream>
#include <vector>
#include <sstream>
#include <string.h>
#include <chrono> 
using namespace std;
using namespace std::chrono;

std::pair<vector<long double>, std::pair<int, int>> getVectorFromFile(string filename);

int main()
{
	std::cout << "Loading matrix from file" << std::endl;

	string file_name = "data_8.txt";
	auto start = high_resolution_clock::now();
	auto input = getVectorFromFile(file_name);
	auto stop = high_resolution_clock::now();
	auto load_file_duration_sec = duration_cast<seconds>(stop - start);
	auto load_file_duration_millis = duration_cast<milliseconds>(stop - start);

	cout << "Load File Duration: " << load_file_duration_sec.count() << " sec" << endl << endl;
	cout << "Load File Duration: " << load_file_duration_millis.count() << " millis" << endl << endl;

	std::cout << "Creating matrix object" << std::endl;
	LA::Matrix mat(input.second.first, input.second.second, &input.first);
	//LA::Matrix mat(5, 5, true);
	std::cout << "Finished Creating matrix object" << std::endl;
	std::cout << "Calculate SVD" << std::endl;
	start = high_resolution_clock::now();
	mat.SVD(mat.GetColNum());

	stop = high_resolution_clock::now();
	std::cout << "Finished Calculating SVD" << std::endl;

	auto SVD_duration_sec = duration_cast<seconds>(stop - start);
	auto SVD_duration_millis = duration_cast<milliseconds>(stop - start);

	cout << "SVD Duration: " << SVD_duration_sec.count() << " sec" << endl << endl;
	cout << "SVD Duration: " << SVD_duration_millis.count() << " millis" << endl << endl;

	auto u = mat.GetU();
	auto s = mat.GetS();
	auto vt = mat.GetVt();

	start = high_resolution_clock::now();
	auto mulRes = u * s * vt;

	auto res = mat.Compare(mulRes, 0.01);
	stop = high_resolution_clock::now();
	auto compare_duration_sec = duration_cast<seconds>(stop - start);
	auto compare_duration_millis = duration_cast<milliseconds>(stop - start);
	std::cout << "Compare result: " << res << std::endl;

	cout << "Compare Duration: " << compare_duration_sec.count() << " sec" << endl << endl;
	cout << "Compare Duration: " << compare_duration_millis.count() << " millis" << endl << endl;

	mat.WriteMatrixToFile("data8", "A");
	mat.GetV().WriteMatrixToFile("data8", "V");
	s.WriteMatrixToFile("data8", "S");
	u.WriteMatrixToFile("data8", "U");
	mulRes.WriteMatrixToFile("data8", "USVt");

}

std::pair<vector<long double>, std::pair<int, int>> getVectorFromFile(string filename)
{
	vector<long double> result;
	ifstream input(filename);
	string s;
	int rows = 0;
	int cols = 0;
	while (getline(input, s))
	{
		cols = 0;
		istringstream iss(s);
		string num;

		while (getline(iss, num, ','))
		{
			result.push_back(stod(num));
			cols++;
		}
		rows++;
	}
	return std::pair<vector<long double>, std::pair<int,int>>(result, pair<int,int>(rows, cols));
}


