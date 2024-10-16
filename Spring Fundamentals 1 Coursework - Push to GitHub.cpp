// matlib24.cpp : This file contains the 'main' function and others to test matrix routines
//

// matlib24.cpp : A simple matrix library program.
// Template Written by Prof Richard Mitchell    6/12/23
// Adapted by: Shyam Popat

#include <string>
#include <iostream>
using namespace std;

const int maxM = 16;	// allows for matrices up to size 4*4

struct myMat {			// structure to store a matrix
	int numRows;		// number of rows
	int numCols;		// number of columns
	int data[maxM];		// data are stored in row order
};

myMat zeroMat(int r, int c) {
	// create a matrix with r rows and c columns, filled with zeros
	myMat m{};			// define matrix
	m.numRows = r;		// set size
	m.numCols = c;
	for (int ct = 0; ct < maxM; ct++) m.data[ct] = 0;	// set elems to 0
	return m;			// return matrix
}

int getIndex(myMat m, int r, int c) {
	// returm index of element in row r, col c of matrix m
	return r * m.numCols + c;
}

myMat matError(string errstr) {
	// if this is called, an error has occured ... output errstr and return 0 size myMat
	cout << "Error : " << errstr << "\n";
	return zeroMat(0, 0);
}

int intError(string errstr) {
	// if this is called, an error has occured ... output errstr and return 0 
	cout << "Error : " << errstr << "\n";
	return 0;
}

myMat mFromStr(string s) {
	// create a matrix from string s
	// string of form "1,2;3,4;5,6"   ; separates rows and , separates columns ... No error check
	int ms;
	if (s.length() > 0) ms = 1; else ms = 0;
	myMat m = zeroMat(ms, ms);						// is s empty create 0*0 matrix, else start with 1*1 matrix
	int mndx = 0;									// used to index into array
	string sub = "";								// sub string - numbers between , or ; set empty
	for (int ct = 0; ct < s.length(); ct++) {		// each char in turn
		if ((s[ct] == ',') || (s[ct] == ';')) {	// if , or ; then sub contains number
			m.data[mndx++] = stoi(sub);				// convert number to integer, put in data array
			sub = "";								// clear sub string
			if (s[ct] == ';') m.numRows++;			// if found ; indicates an extra row
			else if (m.numRows == 1) m.numCols++;	// if , then (if in row 1) increase count of columns
		}
		else sub = sub + s[ct];						// add character to sub string
	}
	if (sub.length() > 0) m.data[mndx++] = stoi(sub);// add last sub string
	return m;
}

void printMat(const char* mess, myMat m) {
	// mess is string to be printed, followed by matrix m
	cout << mess << " = " << "\n";				// print message
	for (int r = 0; r < m.numRows; r++) {		// do each row
		for (int c = 0; c < m.numCols; c++)		// do each column
			cout << m.data[getIndex(m, r, c)] << "\t";	// outputing the element then tab
		cout << "\n";							// output new line at end of row
	}
	cout << "\n";								// and end of Matrix
}


myMat mGetRow(myMat m, int row) {
	// create a matrix from m, having one row
	myMat res = zeroMat(1, m.numCols);		// create a matrix of 1 row
	for (int col = 0; col < m.numCols; col++)		// for each element in row
		res.data[col] = m.data[getIndex(m, row, col)];		// copy col element to res
	return res; //return res
}

myMat mGetCol(myMat m, int col) {
	// create a matrix from m, having one col (done)
	myMat res = zeroMat(m.numRows, 1);	// create a matrix for 1 column
	for (int row = 0; row < m.numRows; row++) // for each element in column
		res.data[row] = m.data[getIndex(m, row, col)]; //copy row element to res
	return res; //return res
}

myMat mSetCol(myMat m, int col, myMat v) {
	// insert v into given col in m
	if (m.numRows != v.numRows) //Checks if the number of rows in matrix m equals the number of rows in matrix v
		return matError("Matrix/Vector should have same number of rows"); //Outputs a message saying the matrix/vector needs to have the same number of rows
	else {
		myMat res = m; //Sets res equal to m
		for (int row = 0; row < m.numRows; row++) //Parses through the number of rows in m (using row as a scope)
			res.data[getIndex(m, row, col)] = v.data[row]; //v row replaces the column in res row
		return res; //returns res
	}
}

int dotProd(myMat v1, myMat v2) {
	int sum = 0; //Sets sum to 0 as an integer
	for (int k = 0; k < v1.numCols; k++) //Parses through the number of columns in v1 (using k as a scope)
		sum = sum + (v1.data[getIndex(v1, 0, k)] * v2.data[getIndex(v2, 0, k)]); //Multiplies matrix v1 and matrix v2 together in the first row and columns increase by 1 (which uses k)
	return sum; //Returns the sum
}	// calculate the dot product of two vectors v1 and v2, each of which could be eitherrow or column vectors

void testVecs(myMat m1, myMat m3) {
	// test vector routines ... get row from m1, col from m3, do dot product of these
	cout << "Testing Vector routines" << "\n";
	printMat("m1 row 0", mGetRow(m1, 0));	// display row 0 of m1
	printMat("m3 col 1", mGetCol(m3, 1));	// display col 1 of m2
	cout << "Dot prod of these is " << dotProd(mGetRow(m1, 0), mGetCol(m3, 1)) << "\n\n";
	cout << "Dot prod of m1 row 1 and m3 row 1 " << dotProd(mGetRow(m1, 0), mGetRow(m3, 1)) << "\n\n";
}

myMat mTranspose(myMat m) {
	myMat res = zeroMat(m.numCols, m.numRows); //Sets up the matrix res
	for (int row = 0; row < m.numRows; row++) //Parses through the number of rows (using m.numRows)
		for (int col = 0; col < m.numCols; col++) //Parses through the number of columns (using m.numCols)
			res.data[getIndex(res, col, row)] = m.data[getIndex(m, row, col)]; //Does the transpose of the matrix
	return res; //res is the matrix which is returned
}

myMat mAdd(myMat m1, myMat m2) {
	myMat res = zeroMat(m1.numRows, m1.numCols); //Sets up the matrix res	
	if (m1.numRows == m2.numRows and m1.numCols == m2.numCols) { //Checks if the number of columns and rows are the same
		for (int row = 0; row < m1.numRows; row++) //Parses through the number of rows (using m.numRows)
			for (int col = 0; col < m1.numCols; col++) //Parses through the number of columns (using m.numCols)
				res.data[getIndex(res, row, col)] = m1.data[getIndex(m1, row, col)] + m2.data[getIndex(m2, row, col)]; //Does the addition of m1 and m2
		return res; //res is the matrix which is returned
	}
	else {
		return matError("\nWrong size!\n"); //Outputs a message if the matrices aren't the same size
	}
}

myMat mScalarMultDiv(myMat m, int s, int isMult) {
	myMat res = zeroMat(m.numRows, m.numCols); //Sets up the matrix res using the number of rows and columns in matrix m
	if (isMult == 1) { //Checks if isMult is equal to 1
		for (int row = 0; row < m.numRows; row++) //Parses through the number of rows (using m.numRows)
			for (int col = 0; col < m.numCols; col++) //Parses through the number of columns (using m.numCols)
				res.data[getIndex(res, row, col)] = m.data[getIndex(m, row, col)] * s; //Multiplies matrix m by s
	}
	else {
		if (s == 0) { //Checks if s is 0 - this is to ensure that the matrix does not divide by 0 to avoid an indefinite loop
			return matError("\nYou cannot divide numbers by 0 in math\n"); //Outputs a message 
		}
		else {
			for (int row = 0; row < m.numRows; row++) //Parses through the number of rows (using m.numRows)
				for (int col = 0; col < m.numCols; col++) //Parses through the number of columns (using m.numCols)
					res.data[getIndex(res, row, col)] = m.data[getIndex(m, row, col)] / s; //Divides matrix m by s
		}
	}
	return res; //Returns matrix res
}

myMat mMult(myMat m1, myMat m2) {
	myMat res = zeroMat(m1.numRows, m2.numCols); //Sets up matrix res using the number of rows in matrix m1 and the number of columns in matrix m2
	if (m1.numCols == m2.numRows) { //Checks if matrix m1 number of rows is not equal to matrix m2 number of columns
		for (int i = 0; i < m1.numRows; i++) { // Iterate through each row of matrix m1
			for (int j = 0; j < m2.numCols; j++) { // Iterate through each column of matrix m2
				myMat rowVector = mGetRow(m1, i); // Get the current row of matrix m1
				myMat colVector = mGetCol(m2, j); // Get the current column of matrix m2
				res.data[getIndex(res, i, j)] = dotProd(rowVector, colVector); // Compute dot product and store in result matrix
			}
		}
		return res; //Returns matrix res
	}
	else {
		return matError("\nWrong sizes!\n"); //Returns an error message saying the matrices are the wrong sizes
	}
}

//Below is for assessment

myMat mSubMat(myMat m, int row, int col) {
	myMat res = zeroMat(m.numRows - 1, m.numCols - 1);
	//Sets up the number of rows and columns for res. We need to take 1 away because we are looking through the indexes of the rows and columns
	int resRow = 0, resCol = 0; // Set resRow and resCol to 0
	for (int r = 0; r < m.numRows; r++) { //Parses through m.numRows using r as a scope
		if (r == row) { //Checks if r is equal to row
			continue;  // Skip the specified row
		}
		for (int c = 0; c < m.numCols; c++) { //Parses through m.numCols using c as a scope
			if (c == col) { //Checks if c is equal to col
				continue;  // Skip the specified column
			}
			res.data[getIndex(res, resRow, resCol)] = m.data[getIndex(m, r, c)]; //Copies the row and the column from matrix m to matrix res
			resCol++; //Moves onto the next column in the result matrix
		}
		// Move to the next row in the result matrix
		// resCol needs to be set to 0 given that a new row has started
		resCol = 0;
		resRow++;
	}
	return res;
	//Returns the matrix res
}

int mDet(myMat m) {
	if (m.numRows != m.numCols) { //Checks if the number of rows is not equal to the columns
		printf("This matrix does not have an equal number of rows and columns!\n"); //A message is output saying that rows is not equal to columns
		return 0; //0 is returned
	}
	else if (m.numRows == 1 and m.numCols == 1) { //Checks if the matrix is 1*1
		return m.data[getIndex(m, 0, 0)]; //Returns the determinant of a 1*1 matrix 
	}
	else if (m.numRows == 2 and m.numCols == 2) { //Checks if the matrix is 2*2
		return m.data[getIndex(m, 0, 0)] * m.data[getIndex(m, 1, 1)] - m.data[getIndex(m, 0, 1)] * m.data[getIndex(m, 1, 0)]; //Returns the determinant of a 2*2 matrix
	}
	else {
		int sum = 0; //Sets integer sum to 0
		for (int col = 0; col < m.numCols; col++) { //Parses through m.numCols using col as a scope
			myMat subMat = mSubMat(m, 0, col); //Equates subMat to the value returned from mSubMat with parameters m, 0 and col
			int sign = (col % 2 == 0) ? 1 : -1; //Swaps the sign around for alternate elements
			sum += sign * m.data[getIndex(m, 0, col)] * mDet(subMat); //Calculates the sum
		}
		return sum; //Returns the determinant of the square matrix bigger than 3*3 
	}
}

myMat mAdj(myMat m) {
	myMat res = zeroMat(m.numRows, m.numCols); //Sets up the number of rows and columns for res
	if (m.numRows != m.numCols) { //Checks if the number of rows is not equal to the columns
		return matError("\nNot a square matrix\n"); //Returns that the matrix needs to have the same number of rows and columns
	}
	else if (m.numRows == 1 and m.numCols == 1) { //Checks if the matrix is 1*1
		res.data[getIndex(res, 0, 0)] = m.data[getIndex(m, 0, 0)]; //gets the value for row 1 col 1 in the 1*1 matrix
		return res; //Returns the 1*1 matrix
	}
	else if (m.numRows == 2 and m.numCols == 2) { //Checks if the matrix is 2*2
		res.data[getIndex(res, 0, 0)] = m.data[getIndex(m, 1, 1)]; //gets the value for row 1 col 1 in the 2*2 matrix (for the cofactors)
		res.data[getIndex(res, 0, 1)] = -m.data[getIndex(m, 1, 0)]; //gets the value for row 1 col 2 in the 2*2 matrix (for the cofactors)
		res.data[getIndex(res, 1, 0)] = -m.data[getIndex(m, 0, 1)]; //gets the value for row 2 col 1 in the 2*2 matrix (for the cofactors)
		res.data[getIndex(res, 1, 1)] = m.data[getIndex(m, 0, 0)]; //gets the value for row 2 col 2 in the 2*2 matrix (for the cofactors)
		res = mTranspose(res); //To get the adjoint, we need to find the transpose of the cofactor matrix
		return res; //Return the 2*2 matrix
	}
	else {
		for (int row = 0; row < m.numRows; row++) { //Parses through m.numRows as row is used as a scope
			for (int col = 0; col < m.numCols; col++) { //Parses through m.numCols as col is used as a scope
				myMat minor = mSubMat(m, row, col); // Calculate the minor (submatrix excluding the current row and column)
				int minorDet = mDet(minor); // Calculate the determinant of the minor
				int cofactor = ((row + col) % 2 == 0) ? minorDet : -minorDet; // Calculate the cofactor (sign changes for alternate elements)
				res.data[getIndex(res, col, row)] = cofactor; // Assign the cofactor to the corresponding position in the adjoint matrix
			}
		}
		return res; //Returns the adjoint of the square matrix bigger than 3*3
	}
}

void testMatEqn(myMat A, myMat b) {
	myMat result_1 = zeroMat(A.numRows, 1); //Sets up 1 row and the number of columns in A for the cramer method (which is result_1)
	myMat result_2 = zeroMat(A.numRows, 1); //Sets up 1 row and the number of columns in A for the adjoint method (which is result_2)
	if (mDet(A) == 0) { //Checks if the dererminant of A is 0
		matError("The determinant of A must be greater than 0!"); //Returns an error message saying the determinant of A must be greater than 0
	}
	else {
		//Cramer method
		for (int row = 0; row < A.numRows; row++) { //Parses through the number of columns using A.numRows
			result_1.data[row] = mDet(mSetCol(A, row, b)) / mDet(A);
			/*
			Store the solution in the result matrix
			Uses the formula mDet(Ai) / mDet(A) but we need to set the column for matrix A and use each row parsed through with matrix b
			So what we have done for mSetCol is that we have replaced the row parsed through in matrix A with matrix b
			We also need to divide by Det(A)
			Note that this finds the value of the individual elements in the matrix
			*/
		}
		printMat("Here are the results from Cramer method: \n", result_1); //Outputs the matrix results
		//Adjoint method
		for (int row = 0; row < A.numRows; row++) { //Parses through the number of columns using A.numRows
			result_2.data[row] = dotProd(mGetRow(mAdj(A), row), b) / mDet(A);
			/*
			Store the solution in the result matrix
			Uses the formula (Adj(A) * b) / Det(A) but we need to use the adjacent of A and get the rows parsed through as well as get the rows for b and use the first column
			We then need to multiply it together using the dot product rule
			We also need to divide by Det(A)
			Note that this finds the value of the individual elements in the matrix
			*/
		}
		printMat("Here are the results from Adjoint method: \n", result_2); //Outputs the matrix results
		//Match verification is below
		for (int row = 0; row < A.numRows; row++) { //Parses through teh number of columns using A.numCols
			if (result_1.data[row] == result_2.data[row]) { //Checks if the values match on each column
				printf("\nBoth values match.\n"); //Outputs a message saying the value matches 
			}
			else {
				matError("\nThe solutions do not match. We do not have the solution unfortunately.\n"); //Outputs a message saying the values don't match
			}
		}
	}
}

void testMatOps(myMat m1, myMat m2, myMat m3) {
	// test matrix operations m1 + m2; m1 + m3 (not poss) m1 + m3'; m1 * m3; m3 * m2
	cout << "Testing Add, Transpose, Multiply routines" << "\n";
	printMat("m1 + m2", mAdd(m1, m2));
	printMat("m1 + m3", mAdd(m1, m3));
	printMat("m1 + m3'", mAdd(m1, mTranspose(m3)));
	printMat("m1 * m3", mMult(m1, m3));
	printMat("m3 * m1", mMult(m3, m1));
	printMat("m1 * m2", mMult(m1, m2));
	printf("\nm1 * [matrix to find] = m2\n");
	testMatEqn(m1, m2);
}

int main() {
	cout << "32010160 Matrix Program\n";						// change to your student number
	myMat m1, m2, m3;											// create matrices
	m1 = mFromStr("8,10;4,10");									// change numbers to those in A from Q1 on web page, as specified on the sheet
	m2 = mFromStr("164;132");									// ditto but B
	m3 = mFromStr("5,4,3,2;2,3,5,5");							// ditto but C	
	printMat("m1", m1);											// display m1
	printMat("m2", m2);											// display m2
	printMat("m3", m3);											// display m3
	testVecs(m1, m3);											// test the vector routines
	testMatOps(m1, m2, m3);										// test the add, transpose and multiply
	return 0;
}

//Note the following:
//m1 = mFromStr("1,2,3;4,-5,6");
//m2 = mFromStr("2,-5,3;7,1,4");
//m3 = mFromStr("2,-5;3,7;1,4");
//These were for the original matrices given
