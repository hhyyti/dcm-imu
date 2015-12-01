# dcm-imu
The DCM-IMU algorithm is designed for fusing low-cost triaxial MEMS gyroscope and accelerometer measurements. An extended Kalman filter is used to estimate attitude in direction cosine matrix (DCM) formation and gyroscope biases online. A variable measurement covariance method is implemented for acceleration measurements to ensure robustness against temporarily non-gravitational accelerations which usually induce errors to attitude estimate in ordinary IMU-algorithms.

If you use the algorithm in any scientific context, please cite: Heikki Hyyti and Arto Visala, “A DCM Based Attitude Estimation Algorithm for Low-Cost MEMS IMUs,” International Journal of Navigation and Observation, vol. 2015, Article ID 503814, 18 pages, 2015. http://dx.doi.org/10.1155/2015/503814

If you would like to use comparison algorithms by Sebastian Madgwick, download them from http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/ and copy c-implementations under c/MahonyAHRS/ (MahonyAHRS.cpp and MahonyAHRS.h) and c/MadgwickAHRS/ (MadgwickAHRS.cpp and MadgwickAHRS.h) folders. The c files have to be renamed as cpp files in order to allow Matlab to compile them correctly. In addition, copy folders @MadgwickAHRS, @MahonyAHRS and quaternion_library into the main folder from the provided Matlab code by Madgwick. These files are not added into this repository as they are provided under GPL licence and this work is under MIT licence.


