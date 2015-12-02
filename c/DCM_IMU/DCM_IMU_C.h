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
 * Default constant values if they are not set in the constructor
 */
#define DEFAULT_g0 9.8189
#define DEFAULT_state {0,0,1,0,0,0}
#define DEFAULT_q_dcm2 (0.1*0.1)
#define DEFAULT_q_gyro_bias2 (0.0001*0.0001)
#define DEFAULT_r_acc2 (0.5*0.5)
#define DEFAULT_r_a2 (10*10)
#define DEFAULT_q_dcm2_init (1*1)
#define DEFAULT_q_gyro_bias2_init (0.1*0.1)

//! DCM_IMU_C class.
/*!
  The DCM-IMU algorithm is designed for fusing low-cost triaxial MEMS gyroscope and accelerometer measurements. An extended Kalman filter is used to estimate attitude in direction cosine matrix (DCM) formation and gyroscope biases online. A variable measurement covariance method is implemented for acceleration measurements to ensure robustness against temporarily non-gravitational accelerations which usually induce errors to attitude estimate in ordinary IMU-algorithms.
  If you use the algorithm in any scientific context, please cite: Heikki Hyyti and Arto Visala, “A DCM Based Attitude Estimation Algorithm for Low-Cost MEMS IMUs,” International Journal of Navigation and Observation, vol. 2015, Article ID 503814, 18 pages, 2015. http://dx.doi.org/10.1155/2015/503814
*/
class DCM_IMU_C {

public:
    //! DCM_IMU_C constructor.
    /*!
      Initializes DCM_IMU_C either with default values or given parameters. All parameters are in SI-units.

      \param Gravity A magnitude of gravity
      \param State An initial state as a array of six doubles, DCM states and bias states.
      \param Covariance A covariance matrix (size of 6x6 doubles, array of 36 doubles in row-major order). If a custom covariance matrix is given, parameters InitialDCMVariance and InitialBiasVariance are not used.
      \param DCMVariance a variance for DCM state update, Q(0,0), Q(1,1), and Q(2,2)
      \param BiasVariance a variance for bias state update, Q(3,3), Q(4,4), and Q(5,5)
      \param InitialDCMVariance an initial variance for DCM state, P(0,0), P(1,1), and P(2,2). If Covariance matrix is given, this parameter is not used.
      \param InitialBiasVariance an initial variance for bias state, P(3,3), P(4,4), and P(5,5). If Covariance matrix is given, this parameter is not used.
      \param MeasurementVariance a constant part of the variance for measurement update, R(0,0), R(1,1), and R(2,2)
      \param MeasurementVarianceVariableGain a gain for the variable part of the variance for measurement update, R(0,0), R(1,1), and R(2,2)
    */
    DCM_IMU_C(const double Gravity = DEFAULT_g0, const double *State = NULL, const double *Covariance = NULL,
        const double DCMVariance = DEFAULT_q_dcm2, const double BiasVariance = DEFAULT_q_gyro_bias2,
        const double InitialDCMVariance = DEFAULT_q_dcm2_init, const double InitialBiasVariance = DEFAULT_q_gyro_bias2_init,
        const double MeasurementVariance = DEFAULT_r_acc2, const double MeasurementVarianceVariableGain = DEFAULT_r_a2);
    
    //! A method to perform update and give new measurements.
    /*!
      This method is used regularly to update new gyroscope and accelerometer measurements into the system. To get best performance of the filter, please calibrate accelerometer and gyroscope readings before sending them into this method. The calibration process is documented in http://dx.doi.org/10.1155/2015/503814
      In addition, please measure the used sample period as accurately as possible for each iteration (delay between current and the last data which was used at the previous update)
      All parameters are in SI-units.

      \param Gyroscope an array of gyroscope measurements (the length is 3 doubles, angular velocities around x, y and z axis).
      \param Accelerometer an array of accelerometer measurements (the length is 3 doubles, accelerations in x, y and z axis).
      \param SamplePeriod A delay between this measurement and the previous measurement in seconds.
    */
    void updateIMU(const double *Gyroscope, const double *Accelerometer, const double SamplePeriod);

    //! A method to query State.
    /*!
      \param State a 6 units length double array where the current state is stored.
    */
    void getState(double *State);

    //! A method to query Covariance.
    /*!
      \param Covariance a 36 units length double array where the current covariance is stored in row-major order.
    */
    void getCovariance(double *Covariance);

    //! A method to query non-gravitational acceleration.
    /*!
      \param a a 3 units length double array where the current non-gravitational acceleration is stored (x, y, and z axis).
    */
    void getNGAcc(double *a);

    //! A method to return the yaw angle.
    /*!
      \return The yaw angle.
    */
    inline double getYaw() { return yaw; }

    //! A method to return the pitch angle.
    /*!
      \return The pitch angle.
    */
    inline double getPitch() { return pitch; }

    //! A method to return the roll angle.
    /*!
      \return The roll angle.
    */
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
