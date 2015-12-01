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

#ifndef DCM_IMU_C_H
#define DCM_IMU_C_H

/**
 * Includes
 */
#include "math.h"
#include <Eigen/Core>
#include <Eigen/Dense>

/**
 * Defines
 */
#define DEFAULT_g0 9.8189
#define DEFAULT_state {0,0,1,0,0,0}
#define DEFAULT_q_dcm2 (0.1*0.1)
#define DEFAULT_q_gyro_bias2 (0.0001*0.0001)
#define DEFAULT_r_acc2 (0.5*0.5)
#define DEFAULT_r_a2 (10*10)
#define DEFAULT_q_dcm2_init (1*1)
#define DEFAULT_q_gyro_bias2_init (0.1*0.1)

class DCM_IMU_C {

public:
    DCM_IMU_C(const double Gravity = DEFAULT_g0, const double *State = NULL, const double *Covariance = NULL,
        const double DCMVariance = DEFAULT_q_dcm2, const double BiasVariance = DEFAULT_q_gyro_bias2,
        const double InitialDCMVariance = DEFAULT_q_dcm2_init, const double InitialBiasVariance = DEFAULT_q_gyro_bias2_init,
        const double MeasurementVariance = DEFAULT_r_acc2, const double MeasurementVarianceVariableGain = DEFAULT_r_a2);
    
    void updateIMU(const double *Gyroscope, const double *Accelerometer, const double SamplePeriod);
    void getState(double *State);
    void getCovariance(double *Covariance);
    void getNGAcc(double *a);
    inline double getYaw() { return yaw; }
    inline double getPitch() { return pitch; }
    inline double getRoll() { return roll; }
    
private:
    double g0;
    Eigen::Matrix<double, 6, 1> x;
    double q_dcm2;
    double q_gyro_bias2;
    double r_acc2;
    double r_a2;
    Eigen::Matrix<double, 3, 1> a;
    double yaw;
    double pitch;
    double roll;
    Eigen::Matrix<double, 6, 6> P;
    Eigen::Matrix<double, 3, 6> H;
    Eigen::Matrix<double, 6, 6> Q;
};


#endif /* DCM_IMU_C_H */
