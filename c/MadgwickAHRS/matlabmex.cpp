//============================================================================
// Copyright (C) 2014, Heikki Hyyti
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//============================================================================

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include "MadgwickAHRS.h"

using namespace std;

extern void _main();

const int numInputArgs  = 2;
const int numOutputArgs = 1;
const int inputCols[2] = {3, 3};

// mexFunction for Matlab
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	if (sizeof(float) * CHAR_BIT != 32) { //Cause error on systems with 64bit float as invSqrt will not work properly!
		mexErrMsgTxt("32 bit floats are required for invSqrt in the included code!");
		return;
	}

	if (sizeof(long) * CHAR_BIT != 32) { //Madgwick's code uses long for 32bit signed integer. It is 64bit in LP64 architecture and the invSqrt miscalculates.
		mexPrintf("Warning: invSqrt does not work with %d bit longs\nPlease change it to int or int32_t in the included code!\n", sizeof(long) * CHAR_BIT);
	}

	// Check to see if we have the correct number of input and output
	// arguments.
	if (nrhs != numInputArgs) {
		mexErrMsgTxt("Incorrect number of input arguments, should be 2");
		return;
	}
	if (nlhs != numOutputArgs) {
		mexErrMsgTxt("Incorrect number of output arguments, should be 1");
		return;
	}

	//use row count of smallest input object
	int rows = mxGetM(prhs[0]);
	for (int i = 0; i < nrhs; ++i) {
		if (rows > mxGetM(prhs[i])) rows = mxGetM(prhs[i]);

		if (inputCols[i] != mxGetN(prhs[i])) {
			mexErrMsgTxt("Incorrect number of input columns, should be 3 and 3");
			return;
		}
	}


	// Create the output structures.
	plhs[0] = mxCreateDoubleMatrix(rows,4,mxREAL); //q

	//inputs for one iteration
	double gyro[3];
	double acc[3];

	//input pointers
	double *gyro_input = mxGetPr(prhs[0]);
	double *acc_input = mxGetPr(prhs[1]);

	//outputs for one iteration
	double q[4];

	//output pointers
	double *q_output = mxGetPr(plhs[0]);

	//(note matrix is organized differently in matlab and c++ ... )
	//input[(rows*collIndex)+rowIndex]

	for (int i = 0; i < rows; ++i) {
		//read input
		for (int j = 0; j < 3; ++j) {
			gyro[j] = gyro_input[(rows*j) + i];
			acc[j] = acc_input[(rows*j) + i];
		}

		//update one iteration
		MadgwickAHRSupdateIMU(gyro[0], gyro[1], gyro[2], acc[0], acc[1], acc[2]);

		//get output from filter
		q[0] = q0;
		q[1] = q1;
		q[2] = q2;
		q[3] = q3;

		//set output
		for (int j = 0; j < 4; ++j) {
			q_output[(rows*j) + i] = q[j];
		}
	}
}

