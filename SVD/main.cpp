
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

	std::cout << "Writing matrices to output file..." << std::endl;
	auto u = mat.GetU();
	auto s = mat.GetS();
	auto vt = mat.GetVt();
	auto mulRes = u * s * vt;

	auto compare_res = mat.Compare(mulRes, 0.01);


	string filename = "data8";
	mat.WriteMatrixToFile(filename, "A");
	mat.GetV().WriteMatrixToFile(filename, "V");
	s.WriteMatrixToFile(filename, "S");
	u.WriteMatrixToFile(filename, "U");
	mulRes.WriteMatrixToFile(filename, "USVt");


	std::ofstream output_file;
	auto res = compare_res == 1 ? "Successful" : "Failed";
	output_file.open(filename.append("_output.txt"), std::ios_base::app);
	output_file << "SVD Results: " << res <<endl << endl;
	output_file << "Performance Results: " << endl;
	output_file << "Load File Duration: " << load_file_duration_sec.count() << " sec" << endl;
	output_file << "Load File Duration: " << load_file_duration_millis.count() << " millis" << endl << endl;
	output_file << endl;
	output_file << "SVD Duration: " << SVD_duration_sec.count() << " sec" << endl;
	output_file << "SVD Duration: " << SVD_duration_millis.count() << " millis" << endl << endl;
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


