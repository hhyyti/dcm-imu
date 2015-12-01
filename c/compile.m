%============================================================================
% Copyright (C) 2015, Heikki Hyyti
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%============================================================================


disp('Compiling DCM-IMU...');
mex -lm -I/usr/include/eigen3 -I./DCM_IMU -I/usr/local/matlab/extern/include ...
    -output ../DCM_IMU_C DCM_IMU/matlabmex.cpp DCM_IMU/DCM_IMU_C.cpp

disp('Compiling Madgwick AHRS...');
mex -lm -I/usr/local/matlab/extern/include -I./MadgwickAHRS ...
    -output ../Madgwick_IMU_C MadgwickAHRS/matlabmex.cpp MadgwickAHRS/MadgwickAHRS.cpp

disp('Compiling Mahony AHRS...');
mex -lm -I/usr/local/matlab/extern/include -I./MahonyAHRS ...
    -output ../Mahony_IMU_C MahonyAHRS/matlabmex.cpp MahonyAHRS/MahonyAHRS.cpp
