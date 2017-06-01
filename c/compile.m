%============================================================================
% Copyright (C) 2017, Heikki Hyyti
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
    -output ../DCM_IMU_C DCM_IMU/matlabmex.cpp DCM_IMU/DCM_IMU_C.cpp;

disp('Compiling DCM-IMU-uc...');
mex -lm -I./DCM_IMU_uc -I/usr/local/matlab/extern/include ...
    -output ../DCM_IMU_uC DCM_IMU_uc/matlabmex.cpp DCM_IMU_uc/DCM_IMU_uC.cpp;

if (exist('MadgwickAHRS/MadgwickAHRS.cpp', 'file') && exist('MadgwickAHRS/MadgwickAHRS.h', 'file'))
    disp('Compiling Madgwick AHRS...');
    mex -lm -I/usr/local/matlab/extern/include -I./MadgwickAHRS ...
        -output ../Madgwick_IMU_C MadgwickAHRS/matlabmex.cpp MadgwickAHRS/MadgwickAHRS.cpp;
else
    disp('"MadgwickAHRS/MadgwickAHRS.cpp" or "MadgwickAHRS/MadgwickAHRS.h" not found.');
    disp('To get reference algorithms, download their implementations from here:');
    disp('http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/');
end

if (exist('MahonyAHRS/MahonyAHRS.cpp', 'file') && exist('MahonyAHRS/MahonyAHRS.h', 'file'))
    disp('Compiling Mahony AHRS...');
    mex -lm -I/usr/local/matlab/extern/include -I./MahonyAHRS ...
        -output ../Mahony_IMU_C MahonyAHRS/matlabmex.cpp MahonyAHRS/MahonyAHRS.cpp;
else
    disp('"MahonyAHRS/MahonyAHRS.cpp" or "MahonyAHRS/MahonyAHRS.h" not found.');
    disp('To get reference algorithms, download their implementations from here:');
    disp('http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/');
end
