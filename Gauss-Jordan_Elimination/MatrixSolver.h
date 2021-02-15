#pragma once
#include "stdafx.h"
#include "MatrixSolver.h"
#include "iostream"
#include <vector>
#include <fstream>
#include <string>
#include <regex>
#include <iomanip>
#include <map>

using namespace std;
using MATRIX = vector<vector<double>>;
using ROW = vector<double>;
using ANSWERS = vector<double>;

class MatrixSolver
{
public:
	MatrixSolver();
	~MatrixSolver();
	ANSWERS answers;
	ANSWERS solveMatrix(string filePath);
	void showAnswers(ANSWERS answers);
	void showExecutionSteps();
	MATRIX apply(string filePath);
};

