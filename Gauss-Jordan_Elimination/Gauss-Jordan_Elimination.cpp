// Gauss-Jordan_Elimination.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "iostream"
#include <vector>
#include <fstream>
#include <string>
#include <regex>
#include <iomanip>
#include <filesystem>
#include "MatrixSolver.h"

using namespace std;
using MATRIX = vector<vector<double>>;
using ROW = vector<double>;

//#define debug

int main()
{
	MatrixSolver matrixSolver = MatrixSolver();
	//matrixSolver.showAnswers(matrixSolver.solveMatrix("matrix.txt"));
	//matrixSolver.eliminate(matrixSolver.readMatrix("matrix.txt"));
	string s = "matrix.txt";
	matrixSolver.apply(s);

//ifdef debug
	//matrixSolver.showExecutionSteps();
//#endif // debug

	int end;
	cin >> end;
    return 0;
}

