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

classdef DCM_IMU < handle
% DCM_IMU Implementation of Hyyti's IMU algorithm
%
%   If you use the algorithm in any scientific context, please cite: 
%   Heikki Hyyti and Arto Visala, “A DCM Based Attitude Estimation Algorithm for Low-Cost MEMS IMUs,” 
%   International Journal of Navigation and Observation, vol. 2015, Article ID 503814, 18 pages, 2015. 
%   http://dx.doi.org/10.1155/2015/503814  
%
%   Date          Author          Notes
%   1/12/2015     Heikki Hyyti    Initial release

    %% Public properties
    properties (Access = public)
        g0 = 9.8189;                % gravitation around Helsinki, Finland (change according to your area)
        state = [0 0 1 0 0 0]';     % States are lowest row of rotation matrix and gyroscope x y and z biases
                                    % (C_31, C_32, C_33, w_b1, w_b2, w_b3)
        q_dcm2 = 0.1^2;             % estimated variance of dcm states (gyro variance per second)
        q_gyro_bias2 = 0.0001^2;    % very small number to make bias change slowly
        r_acc2 = 0.5^2;             % variance of calibrated accelerometer (g-component)
        r_a2 = 10^2;                % large variance for some unknown acceleration (acc = a + g)
        q_dcm2_init = 1^2;          % initial variance of dcm states (for attitude estimation)
        q_gyro_bias2_init = 0.1^2;  % initial variance of bias states (for bias estimator)
        a = zeros(3,1);             % estimated non-gravitational accelerations
        yaw = 0;                    % Yaw angle around z axis (in ZYX convention)
        pitch = 0;                  % Pitch angle around y axis
        roll = 0;                   % Roll angle around x axis
        P = [];                     % estimate covariance (these are initialized in constructor below)
        H = [];                     % observation model (static)
        Q = [];                     % proces noise covariance (static part)
    end

    %% Public methods
    methods (Access = public)
        function obj = DCM_IMU(varargin)
            updateP = true;
            for i = 1:2:nargin
                if  strcmp(varargin{i}, 'Gravity'), obj.g0 = varargin{i+1};
                elseif  strcmp(varargin{i}, 'State'), obj.state = varargin{i+1};
                elseif  strcmp(varargin{i}, 'Covariance'), obj.P = varargin{i+1}; updateP = false;
                elseif  strcmp(varargin{i}, 'DCMVariance'), obj.q_dcm2 = varargin{i+1};
                elseif  strcmp(varargin{i}, 'BiasVariance'), obj.q_gyro_bias2 = varargin{i+1};
                elseif  strcmp(varargin{i}, 'InitialDCMVariance'), obj.q_dcm2_init = varargin{i+1};
                elseif  strcmp(varargin{i}, 'InitialBiasVariance'), obj.q_gyro_bias2_init = varargin{i+1};                    
                elseif  strcmp(varargin{i}, 'MeasurementVariance'), obj.r_acc2 = varargin{i+1};
                elseif  strcmp(varargin{i}, 'MeasurementVarianceVariableGain'), obj.r_a2 = varargin{i+1};
                else error('Invalid argument');
                end
            end;
            
            if (updateP), obj.P = [obj.q_dcm2_init*eye(3), zeros(3,3); zeros(3,3), obj.q_gyro_bias2_init*eye(3)]; end;
            obj.H = [eye(3)*obj.g0, zeros(3,3)];
            obj.Q = [obj.q_dcm2*eye(3), zeros(3,3); zeros(3,3) obj.q_gyro_bias2*eye(3)];
        end
        function obj = UpdateIMU(obj, Gyroscope, Accelerometer, SamplePeriod)
            x = obj.state;
            Q_ = SamplePeriod^2 * obj.Q; %Process noise covariance with time dependent noise
            
            % control input (angular velocities from gyroscopes)
            if (size(Gyroscope,1) == 3), u = Gyroscope;
            else u = Gyroscope';
            end
            
            % "rotation operators"
            C3X = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
            UX = [0 -(u(3)-x(6)) u(2)-x(5); 
                u(3)-x(6) 0 -(u(1)-x(4)); 
                -(u(2)-x(5)) u(1)-x(4) 0];

            % Model generation
            A = [zeros(3,3) -SamplePeriod*C3X; zeros(3,6)];
            B = [SamplePeriod*C3X; zeros(3,3)];
            F = eye(6) + [-SamplePeriod*UX, -SamplePeriod*C3X; zeros(3,6)];

            % Kalman a priori prediction
            x_predict = x + A*x + B*u;
            P_predict = F * obj.P * F' + Q_;

            % measurements/observations (acceleromeres)
            if (size(Accelerometer,1) == 3), z = Accelerometer;
            else z = Accelerometer';
            end            
            
            % recompute R using the error between acceleration and the model of g 
            % (estimate of the magnitude of a0 in a = a0 + g)
            a_predict = z - x_predict(1:3)*obj.g0;
            a_len = sqrt(a_predict'*a_predict);
            R = (a_len*obj.r_a2 + obj.r_acc2)*eye(3);

            % Kalman innovation
            y = z - obj.H*x_predict;
            S = obj.H * P_predict * obj.H' + R;

            % Kalman gain
            K = P_predict * obj.H' / S;

            % update a posteriori
            x = x_predict + K * y;

            % update a posteriori covariance
            IKH = eye(6) - K*obj.H;
            obj.P = IKH * P_predict * IKH' + K * R * K'; % for using any K

            % normalization of x & P (divide by DCM vector length)
            dcm_vector_length = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
            J_33 = [x(2)^2 + x(3)^2,    -x(1)*x(2),         -x(1)*x(3); ...
                    -x(1)*x(2),         x(1)^2 + x(3)^2,    -x(2)*x(3); ...
                    -x(1)*x(3),         -x(2)*x(3),         x(1)^2 + x(2)^2];        
            J = [ J_33 / (dcm_vector_length^3), zeros(3,3); zeros(3,3), eye(3)];

            % Laplace approximation of normalization function for x to P, J = Jacobian(f,x)
            % P_new = E[J*(x-x0)*(x-x0)'*J'] = J*E[(x-x0)*(x-x0)']*J' = J*P*J'
            obj.P = J*obj.P*J';
            x(1:3) = x(1:3) ./ dcm_vector_length;
            obj.state = x;

            
            % compute Euler angles (not exactly a part of the extended Kalman filter)
            % yaw integration through full rotation matrix
            u_nb = u - x(4:6);

            cy = cos(obj.yaw); %old angles (last state before integration)
            sy = sin(obj.yaw);
            cp = cos(obj.pitch);
            sp = sin(obj.pitch);
            cr = cos(obj.roll);
            sr = sin(obj.roll);

            % compute needed parts of rotation matrix R
            R11 = cy*cp;
            R12 = cy*sp*sr-sy*cr;
            R13 = cy*sp*cr+sy*sr;
            R21 = sy*cp;
            R22 = sy*sp*sr+cy*cr;
            R23 = sy*sp*cr-cy*sr;
            
            % update needed parts of R for yaw computation
            R11_new = R11 + SamplePeriod*(u_nb(3)*R12 - u_nb(2)*R13);
            R21_new = R21 + SamplePeriod*(u_nb(3)*R22 - u_nb(2)*R23);
            obj.yaw = atan2(R21_new,R11_new);

            %compute new pitch and roll angles from a posteriori states
            obj.pitch = asin(-x(1));
            obj.roll = atan2(x(2),x(3));

            % save the estimated non-gravitational acceleration
            obj.a = z - x(1:3)*obj.g0; % acceleration estimate (g reduced)    
        end
    end
end
