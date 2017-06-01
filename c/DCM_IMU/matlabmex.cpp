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
#include "DCM_IMU_C.h"

using namespace std;

extern void _main();

const int numInputArgs  = 3;
const int numOutputArgs = 4;
const int inputCols[3] = {3, 3, 1};

// mexFunction for Matlab
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	// Check to see if we have the correct number of input and output
	// arguments.
	if (nrhs != numInputArgs) {
		mexErrMsgTxt("Incorrect number of input arguments, should be 3");
		return;
	}
	if (nlhs != numOutputArgs) {
		mexErrMsgTxt("Incorrect number of output arguments, should be 4");
		return;
	}

	//use row count of smallest input object
	int rows = mxGetM(prhs[0]);
	for (int i = 0; i < nrhs; ++i) {
		if (rows > mxGetM(prhs[i])) rows = mxGetM(prhs[i]);

		if (inputCols[i] != mxGetN(prhs[i])) {
			mexErrMsgTxt("Incorrect number of input columns, should be 3, 3 and 1");
			return;
		}
	}

	// Create DCM_IMU with default parameters
	DCM_IMU_C imu;

	// Create the output structures.
	plhs[0] = mxCreateDoubleMatrix(rows,6,mxREAL); //x
	plhs[1] = mxCreateDoubleMatrix(rows,3,mxREAL); //ypr
	plhs[2] = mxCreateDoubleMatrix(rows,3,mxREAL); //a
	plhs[3] = mxCreateDoubleMatrix(rows,6,mxREAL); //diag(P)

	//inputs for one iteration
	double gyro[3];
	double acc[3];
	double dt;

	//input pointers
	double *gyro_input = mxGetPr(prhs[0]);
	double *acc_input = mxGetPr(prhs[1]);
	double *time_input = mxGetPr(prhs[2]);

	//outputs for one iteration
	double state[6];
	double ypr[3];
	double a[3];
	double P[36];


	//output pointers
	double *state_output = mxGetPr(plhs[0]);
	double *ypr_output = mxGetPr(plhs[1]);
	double *a_output = mxGetPr(plhs[2]);
	double *diagP_output = mxGetPr(plhs[3]);


	//(note matrix is organized differently in matlab and c++ ... )
	//input[(rows*collIndex)+rowIndex]

	for (int i = 0; i < rows; ++i) {
		//read input
		for (int j = 0; j < 3; ++j) {
			gyro[j] = gyro_input[(rows*j) + i];
			acc[j] = acc_input[(rows*j) + i];
		}
		dt = time_input[i];

		//update one iteration
		imu.updateIMU(gyro, acc, dt);

		//get output from filter
		imu.getState(state);
		imu.getCovariance(P);
		imu.getNGAcc(a);
		ypr[0] = imu.getYaw();
		ypr[1] = imu.getPitch();
		ypr[2] = imu.getRoll();

		//set output
		for (int j = 0; j < 6; ++j) {
			state_output[(rows*j) + i] = state[j];
			diagP_output[(rows*j) + i] = P[7*j];
		}
		for (int j = 0; j < 3; ++j) {
			ypr_output[(rows*j) + i] = ypr[j];
			a_output[(rows*j) + i] = a[j];
		}

	}
}

