#include "stdafx.h"
#include "MatrixSolver.h"
#include "iostream"
#include <vector>
#include <fstream>
#include <string>
#include <regex>
#include <iomanip>
#include <map>

MATRIX initialMatrix;
MATRIX echelonFormMatrix;
ANSWERS answers;
bool solved = false;

MatrixSolver::MatrixSolver()
{
}


MatrixSolver::~MatrixSolver()
{
}

#define debug

int getPivotIndex(ROW r) {
	for (int col = 0; col < r.size(); col++) {
		if (r[col] != 0)
			return col;
	}
	return r.size();
}

MATRIX sortMatrix(MATRIX mat) {
	MATRIX out = mat;
	bool run = true;
	while (run) {
		int swapCounter = 0;
		for (int row = 0; row < out.size(); row++) {
			if (row < out.size() - 1) {
				if (getPivotIndex(out[row]) > getPivotIndex(out[row + 1])) {
					ROW thisRow = out[row];
					out[row] = out[row + 1];
					out[row + 1] = thisRow;
					swapCounter++;
				}
			}
		}
		if (swapCounter == 0)
			run = false;
	}
	return out;
}

// Reads matrix data from file
MATRIX getMatrix(string filePath) {
	vector<vector<double>> out;
	ifstream infile(filePath);
	string line;
	while (getline(infile, line))
	{
		//cout << line << endl;
		regex rx(R"((?:^|\s)([+-]?[[:digit:]]+(?:\.[[:digit:]]+)?)(?=$|\s))");
		smatch m;
		string str = line;
		vector<double> row;
		while (regex_search(str, m, rx))
		{
			// string to int
			double d = stod(m[1]);
			row.push_back(d);
			// To next
			str = m.suffix().str();
		}
		out.push_back(row);
	}
	initialMatrix = out;
	return out;
}

void printMatrix(MATRIX mat, int space = 0) {
	for (int s = 0; s < space; s++)
		cout << "	";
	cout << "Matrix: " << endl;
	for (int row = 0; row < mat.size(); row++) {
		for(int s = 0; s < space; s++)
			cout << "	";
		for (int col = 0; col < mat[row].size(); col++)
		{
			cout << mat[row][col] << " ";
		}
		cout << endl;
	}
}

// returns value from top
double getTopValue(MATRIX *mat, int row, int col) {
	return (*mat)[row - 1][col];
}

ROW getTopROW(MATRIX *mat, int row) {
	return (*mat)[row - 1];
}

double getAt(MATRIX *mat, int row, int col) {
	return (*mat)[row][col];
}

void addToRow(MATRIX *mat, int row, ROW r) {
	for (int col = 0; col < (*mat)[row].size(); col++) {
		(*mat)[row][col] += r[col];
	}
}

void multiplyROW(ROW *r, double value) {
	for (int col = 0; col < (*r).size(); col++)
		(*r)[col] *= value;
}

bool rowEchelonForm(MATRIX mat) {
#ifdef debug
	cout << "rowEchelonForm : " << endl;
#endif // debug

	for (int col = 0; col < mat[0].size() - 1; col++) // mat[0].size() - 1 to ignore the last column (after =)
	{
		// Counter of zeros in this column
		for (int row = mat.size() - 1; row >= col + 1; row--)
		{
			if (mat[row][col] != 0)
				return false;
#ifdef debug
			// counts rows from the bottom of the matrix
			int counter = (row - (mat.size() - 1)) *(-1);
			cout << mat[row][col] << ", " << counter << ",, ";
#endif // debug
		}
	}
	return true;
}

// Use Gauss-Jordan Elimination to obtain Row Echelon Form
MATRIX eliminate(MATRIX mat) {
	MATRIX out = mat;
	for (int col = 0; col < mat.size() - 1; col++)
	{
#ifdef debug
		cout << "col = " << col << endl;
#endif // debug

		ROW initialTopRow = out[col];
		double topCell = getTopValue(&out, col + 1, col);
		//int localCol = col;
		for (int row = col + 1; row < mat.size(); row++)
		{
			if (rowEchelonForm(out))
				return out;
			if (topCell != 0) {
#ifdef debug
				cout << "	row = " << row << endl;
#endif // debug
				// Add sth to a row to obtain 0
				//ROW topRow = getTopROW(&out, row);
				ROW topRow = initialTopRow;
				double thisCell = getAt(&out, row, col);
				//double topCell = getTopValue(&out, row, col);
				double multiplyBy = (-thisCell) / topCell;

#ifdef debug
				cout << "		topRow = ";
				for (int r = 0; r < topRow.size(); r++)
					cout << topRow[r] << ", ";
				cout << endl;

				cout << "		thisCell = " << thisCell << endl;
				cout << "		topCell = " << topCell << endl;
				cout << "		multiplyBy = " << multiplyBy << endl;
				printMatrix(out);
#endif // debug
				multiplyROW(&topRow, multiplyBy);
				addToRow(&out, row, topRow);
			}
		}
	}
	if (!rowEchelonForm(out)) {
		cout << "After applying Gauss-Jordan Elimination matrix if not in a row echelon form!" << endl;
		printMatrix(out);
	}
	echelonFormMatrix = out;
	return out;
}


void MatrixSolver::showAnswers(ANSWERS answers) {
	cout << endl << "	Answers: " << endl;
	for (int i = 0; i < answers.size(); i++)
	{
		cout << "X" << i << " = " << answers[i] << endl;
	}
}

// Get unknown variables if eliminatedMatrix is in echelon form
ANSWERS getAnswers(MATRIX eliminatedMatrix) {

#ifdef debug
	printMatrix(eliminatedMatrix);
#endif // debug

	// Set length = number of columns - 1 (do not count values after equals) 
	ANSWERS out = ANSWERS(eliminatedMatrix[0].size() - 1);
	for (int row = eliminatedMatrix.size() - 1; row >= 0; row--)
	{
			double afterTheUnknown = 0; // 0 + 0 + 2x (+ 5*known) = 10
			// Get sum of numbers after the unknown and before equals 
			// ( col == row to proceed with a triangle)
			for (int i = row + 1; i < eliminatedMatrix[0].size() - 1; i++)
				afterTheUnknown += eliminatedMatrix[row][i] * out[i];
			double nearTheUnknown = eliminatedMatrix[row][row];
			// Value after equals in this row
			double afterEquals = eliminatedMatrix[row][eliminatedMatrix[0].size() - 1];
			if (nearTheUnknown != 0) {
				out[row] = (afterEquals - afterTheUnknown) / nearTheUnknown; // 2x = 8 -> x  8 / 2

#ifdef debug
				cout << "afterEquals = " << afterEquals << endl;
				cout << "afterTheUnknown =  " << afterTheUnknown << endl;
				cout << "nearTheUnknown = " << nearTheUnknown << endl;
				cout << endl << "thisAnswer = " << out[row] << endl;

#endif // debug
			}
	}

	answers = out;
	return out;
}

// Solve using Gauss-Jordan Elimination, returns ROW of answers
ANSWERS MatrixSolver::solveMatrix(string filePath) {
	ANSWERS out = getAnswers(eliminate(getMatrix(filePath)));

#ifdef debug
	cout << endl << "	Answers : " << endl;
	for (int i = 0; i < out.size(); i++)
		cout << "	" << out[i] << endl;
#endif // debug
	solved = true;
	return out;
}

void MatrixSolver::showExecutionSteps() {
	// Show execution steps
	if(solved)
	{
		cout << endl;
		cout << endl << "	Execution Steps : " << endl;
		cout << endl << "		Initial Matrix : " << endl;
		printMatrix(initialMatrix, 3);
		cout << endl << "			Echelon Form Matrix (Eliminated Matrix) : " << endl;
		printMatrix(echelonFormMatrix, 4);
	}
}

// eliminate Reduced Row Echelon Form
MATRIX eliminateRREF(MATRIX mat) {
	MATRIX out = mat;
	// 1. Set pivots equal to (1)
	for (int row = 0; row < out.size(); row++) {
		ROW newRow = out[row];
		int index = getPivotIndex(newRow);
		if (index < newRow.size()) {
			double pivotValue = newRow[index];
			multiplyROW(&newRow, 1 / pivotValue);
			out[row] = newRow;
		}
	}

	// 2. Make zeros on top of each pivot
	for (int row = 1; row < out.size(); row++) {
		int index = getPivotIndex(out[row]);
		if (index < out[0].size())
			for (int topRow = row - 1; topRow >= 0; topRow--) {
				ROW thisRow = out[row];
				ROW nextRow = out[topRow];
				double thisValue = thisRow[index];
				if (thisValue != 0) {
					double nextValue = nextRow[index];
					multiplyROW(&thisRow, -(nextValue / thisValue));
					addToRow(&out, topRow, thisRow);
				}
			}
	}

	return out;
}

MATRIX MatrixSolver::apply(string filePath) {
	MATRIX out;
	cout << endl << "	Output: " << endl << endl;
	MATRIX eliminatedMatrix = sortMatrix(eliminate(getMatrix(filePath)));
	printMatrix(eliminatedMatrix);
	out = sortMatrix(eliminateRREF(eliminatedMatrix));
	cout << endl << "	Final Output (Matrix in RREF): " << endl;
	printMatrix(out, 2);
	return out;
}


