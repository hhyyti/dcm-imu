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

disp(' ');
disp('If you use the algorithm in any scientific context, please cite:');
disp('Heikki Hyyti and Arto Visala, "A DCM Based Attitude Estimation Algorithm');
disp('for Low-Cost MEMS IMUs," International Journal of Navigation and Observation,');
disp('vol. 2015, Article ID 503814, 18 pages, 2015. http://dx.doi.org/10.1155/2015/503814');
disp(' ');

clear all;
close all;
pause(0.1);

%% Configuration

% use c/compile.m before setting this to true (tested only with linux)
use_C_versions = false; 

% fetch GPL licensed algortihms in matlab and C from here before use: 
% http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/ 
% Add the C code to c folder and matlab functions in this folder 
% (please add also the quaternion library used by Madgwick).
% I apologize this mess. The Madgwick's GPS licenced code is separated from 
% this MIT licenced product to avoid infecting this code with GPL licence.
use_comparisonAlgorithms = false; 

% Add bias to calibrated gyro values (same for each axis)
addedGyroBias = 0; 

% load measurement data
data = load('allData.mat');

if (~exist('figures', 'dir')); mkdir('figures'); end;

ILtime = data.InertiaLink.tv_sec + (data.InertiaLink.tv_msec/1000);
SFtime = data.SparkFun6DOF.tv_sec + (data.SparkFun6DOF.tv_usec / 1000000);
Ktime = data.Kuka.tv_sec + (data.Kuka.tv_usec/1000000);

firstCommonTime = max([ILtime(1), SFtime(1), Ktime(1)]);
lastCommonTime = min([ILtime(end), SFtime(end), Ktime(end)]);
commonDuration = lastCommonTime - firstCommonTime;

ILtime = ILtime - firstCommonTime;
SFtime = SFtime - firstCommonTime;
Ktime = Ktime - firstCommonTime;

ILdeltaTime = [0; ILtime(2:end) - ILtime(1:end-1)];
SFdeltaTime = [0; SFtime(2:end) - SFtime(1:end-1)];
KdeltaTime = [0; Ktime(2:end) - Ktime(1:end-1)];

%start yaw from zero at time 468s (start of tests after calibration sequences in example data set)
StartIdx = find(ILtime > 464, 1, 'first');
StopIdx = find(ILtime > 774, 1, 'first');
time = ILtime(StartIdx:StopIdx);

%getting reference trajectory
R11 = data.Kuka.data_msrCartPos01;
R12 = data.Kuka.data_msrCartPos02;
R13 = data.Kuka.data_msrCartPos03;
x = data.Kuka.data_msrCartPos04;
R21 = data.Kuka.data_msrCartPos05;
R22 = data.Kuka.data_msrCartPos06;
R23 = data.Kuka.data_msrCartPos07;
y = data.Kuka.data_msrCartPos08;
R31 = data.Kuka.data_msrCartPos09;
R32 = data.Kuka.data_msrCartPos10;
R33 = data.Kuka.data_msrCartPos11;
z = data.Kuka.data_msrCartPos12;

[yaw, pitch, roll] = yawpitchroll(R11, R21, R31, R32, R33);
yaw = un_modulo(yaw, 2*pi); %remove 2pi jumps and make yaw continuous
yaw = yaw - yaw(1); %starting yaw from zero (same as IMU estimates)

vel_x = [0; x(2:end)-x(1:end-1)] ./ KdeltaTime;
vel_y = [0; y(2:end)-y(1:end-1)] ./ KdeltaTime;
vel_z = [0; z(2:end)-z(1:end-1)] ./ KdeltaTime;

% compute reference angular velocities using the KUKA robot 
% (derivate between positions)
gyro = zeros(length(R11),3);
dt = median(Ktime(2:end)-Ktime(1:end-1)); % the most frequently used discretation time
for i = 2:length(R11)
    R_last = [R11(i-1), R12(i-1), R13(i-1); ...
              R21(i-1), R22(i-1), R23(i-1); ...
              R31(i-1), R32(i-1), R33(i-1)];
    R_cur = [R11(i), R12(i), R13(i); ...
             R21(i), R22(i), R23(i); ...
             R31(i), R32(i), R33(i)];
         
    Wh = (R_last'*R_cur - eye(3));
    W = 1/dt * Wh;
    gyro(i,1) = (W(3,2) - W(2,3))/2;
    gyro(i,2) = (W(1,3) - W(3,1))/2;
    gyro(i,3) = (W(2,1) - W(1,2))/2;
end

% Resample reference data of KUKA robot to InertiaLink time
warning('off','interpolation:interpolation:noextrap');
yaw_in_IL = resample(timeseries(yaw, Ktime),time);
yaw_in_IL = yaw_in_IL.Data;
yaw_in_IL = yaw_in_IL - yaw_in_IL(1); %starting yaw from zero (same as IMU estimates)
pitch_in_IL = resample(timeseries(pitch, Ktime),time);
pitch_in_IL = pitch_in_IL.Data;
roll_in_IL = resample(timeseries(roll, Ktime),time);
roll_in_IL = roll_in_IL.Data;

%computing ypr for InertiaLink internal estimates (Rotation matrix R = -M')
[yaw_IL, pitch_IL, roll_IL] = yawpitchroll(-data.InertiaLink.M11, ...
    data.InertiaLink.M12, -data.InertiaLink.M13, -data.InertiaLink.M23, ...
    -data.InertiaLink.M33);
yaw_IL = un_modulo(yaw_IL, 2*pi); % remove 2pi jumps and make yaw continuous
yaw_IL = yaw_IL - yaw_IL(1); % starting yaw from zero (same as IMU estimates)

%% calibrate IMU data
imusWithKukaCalibration;

acc1 = [acc1_x, acc1_y, acc1_z];
acc2 = [acc2_x, acc2_y, acc2_z];

gyro1 = [w1_x, w1_y, w1_z];
gyro2 = [w2_x, w2_y, w2_z];

% add constant bias for testing purposes (if addedGyroBias is set above)
gyro1 = gyro1 + addedGyroBias*ones(size(gyro1));
gyro2 = gyro2 + addedGyroBias*ones(size(gyro2));

%% DCM IMU results for both imu datas
tic();
if (use_C_versions) 
    [x_hist1, ypr_hist1, a_hist1, P_diag_hist1] = ...
        DCM_IMU_C(gyro1, acc1, SFdeltaTime);
else
    IMU_DCM1 = DCM_IMU();
    x_hist1 = zeros(length(SFtime), 6);
    ypr_hist1 = zeros(length(SFtime), 3);
    P_diag_hist1 = zeros(length(SFtime), 6);
    for t = 1:length(SFtime)
        IMU_DCM1.UpdateIMU(gyro1(t,:), acc1(t,:), SFdeltaTime(t));	% gyroscope units must be radians
        x_hist1(t, :) = IMU_DCM1.state';
        ypr_hist1(t, :) = [IMU_DCM1.yaw, IMU_DCM1.pitch, IMU_DCM1.roll];
        P_diag_hist1(t, :) = diag(IMU_DCM1.P)';
    end
end
tDCM1 = toc()

tic();
if (use_C_versions) 
    [x_hist2, ypr_hist2, a_hist2, P_diag_hist2] = ...
        DCM_IMU_C(gyro2, acc2, ILdeltaTime);
else
    IMU_DCM2 = DCM_IMU();
    x_hist2 = zeros(length(ILtime), 6);
    ypr_hist2 = zeros(length(ILtime), 3);
    P_diag_hist2 = zeros(length(ILtime), 6);
    for t = 1:length(ILtime)
        IMU_DCM2.UpdateIMU(gyro2(t,:), acc2(t,:), ILdeltaTime(t));	% gyroscope units must be radians
        x_hist2(t, :) = IMU_DCM2.state';
        ypr_hist2(t, :) = [IMU_DCM2.yaw, IMU_DCM2.pitch, IMU_DCM2.roll];
        P_diag_hist2(t, :) = diag(IMU_DCM2.P)';
    end
end
tDCM2 = toc()

%% computing Madgwick and Mahony imu algorithms for comparison
if (use_comparisonAlgorithms)
    addpath('quaternion_library');      % include quaternion library

    tic();
    if (use_C_versions) 
        quat_mad_B = Madgwick_IMU_C(gyro1, acc1);
    else
        IMU_mad_B = MadgwickAHRS('SamplePeriod', 1/150, 'Beta', 0.1);
        quat_mad_B = zeros(length(SFtime), 4);
        for t = 1:length(quat_mad_B)
            IMU_mad_B.UpdateIMU(gyro1(t,:), acc1(t,:));	% gyroscope units must be radians
            quat_mad_B(t, :) = IMU_mad_B.Quaternion;
        end
    end
    tIMU_mad_B = toc()

    tic();
    if (use_C_versions) 
        quat_mah_B = Mahony_IMU_C(gyro1, acc1);
    else
        IMU_mah_B = MahonyAHRS('SamplePeriod', 1/150, 'Kp', 0.5);
        quat_mah_B = zeros(length(SFtime), 4);
        for t = 1:length(quat_mah_B)
            IMU_mah_B.UpdateIMU(gyro1(t,:), acc1(t,:));	% gyroscope units must be radians
            quat_mah_B(t, :) = IMU_mah_B.Quaternion;
        end
    end

    tIMU_mah_B = toc()

    tic();
    if (use_C_versions) 
        quat_mad_A = Madgwick_IMU_C(gyro2, acc2);
    else
        IMU_mad_A = MadgwickAHRS('SamplePeriod', 1/150, 'Beta', 0.1);
        quat_mad_A = zeros(length(ILtime), 4);
        for t = 1:length(quat_mad_A)
            IMU_mad_A.UpdateIMU(gyro2(t,:), acc2(t,:));	% gyroscope units must be radians
            quat_mad_A(t, :) = IMU_mad_A.Quaternion;
        end
    end
    tIMU_mad_A = toc()

    tic();
    if (use_C_versions) 
        quat_mah_A = Mahony_IMU_C(gyro2, acc2);
    else
        IMU_mah_A = MahonyAHRS('SamplePeriod', 1/150, 'Kp', 0.5);
        quat_mah_A = zeros(length(ILtime), 4);
        for t = 1:length(quat_mah_A)
            IMU_mah_A.UpdateIMU(gyro2(t,:), acc2(t,:));	% gyroscope units must be radians
            quat_mah_A(t, :) = IMU_mah_A.Quaternion;
        end
    end
    tIMU_mah_A = toc()

    % Plot algorithm output as Euler angles
    % The first and third Euler angles in the sequence (phi and psi) become
    % unreliable when the middle angles of the sequence (theta) approaches ~90
    % degrees. This problem commonly referred to as Gimbal Lock.
    % See: http://en.wikipedia.org/wiki/Gimbal_lock
    
    % use conjugate for sensor frame relative to Earth and convert to degrees.
    euler_mad_B = quatern2euler(quaternConj(quat_mad_B)) * (180/pi);	
    euler_mah_B = quatern2euler(quaternConj(quat_mah_B)) * (180/pi);	
    euler_mad_A = quatern2euler(quaternConj(quat_mad_A)) * (180/pi);	
    euler_mah_A = quatern2euler(quaternConj(quat_mah_A)) * (180/pi);	
else
    euler_mad_B = nan(length(gyro1),3);
    euler_mah_B = nan(length(gyro1),3);
    euler_mad_A = nan(length(gyro2),3);
    euler_mah_A = nan(length(gyro2),3);
    
    disp('To get reference algorithms, download their implementations from here:');
    disp('http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/');
end


%% resampling & remove 2pi jumps
warning('off','interpolation:interpolation:noextrap');

%resample SparkFun data to InertiaLink time
yaw_sf_in_IL = resample(timeseries(un_modulo(ypr_hist1(:,1),2*pi), SFtime),time);
yaw_sf_in_IL = yaw_sf_in_IL.Data*180/pi;
yaw_sf_in_IL = yaw_sf_in_IL - yaw_sf_in_IL(1); %starting yaw from zero (same as IMU estimates)
pitch_sf_in_IL = resample(timeseries(ypr_hist1(:,2), SFtime),time);
pitch_sf_in_IL = pitch_sf_in_IL.Data*180/pi;
roll_sf_in_IL = resample(timeseries(ypr_hist1(:,3), SFtime),time);
roll_sf_in_IL = roll_sf_in_IL.Data*180/pi;


% start yaw from zero at the beginning of test sequence
yaw_IL = yaw_IL - yaw_IL(StartIdx);

%yaw pitch roll errors (only InertiaLink data)
ypr_hist2(:,1) = un_modulo(ypr_hist2(:,1), 2*pi); %remove 2pi jumps and make continuous
ypr_hist2(:,1) = ypr_hist2(:,1) - ypr_hist2(StartIdx,1); %DCM

euler_mad_B(:,3) = un_modulo(euler_mad_B(:,3), 360); %remove 2pi jumps and make continuous
r_data = resample(timeseries(euler_mad_B(:,1), SFtime),ILtime);
p_data = resample(timeseries(euler_mad_B(:,2), SFtime),ILtime);
y_data = resample(timeseries(euler_mad_B(:,3), SFtime),ILtime);
euler_mad_B = [r_data.data, p_data.data, y_data.data];
euler_mad_B(:,3) = euler_mad_B(:,3) - euler_mad_B(StartIdx,3); %Madgwick

euler_mah_B(:,3) = un_modulo(euler_mah_B(:,3), 360); %remove 2pi jumps and make continuous
r_data = resample(timeseries(euler_mah_B(:,1), SFtime),ILtime);
p_data = resample(timeseries(euler_mah_B(:,2), SFtime),ILtime);
y_data = resample(timeseries(euler_mah_B(:,3), SFtime),ILtime);
euler_mah_B = [r_data.data, p_data.data, y_data.data];
euler_mah_B(:,3) = euler_mah_B(:,3) - euler_mah_B(StartIdx,3); %Mahony

euler_mad_A(:,3) = un_modulo(euler_mad_A(:,3), 360); %remove 2pi jumps and make continuous
euler_mad_A(:,3) = euler_mad_A(:,3) - euler_mad_A(StartIdx,3); %Madgwick

euler_mah_A(:,3) = un_modulo(euler_mah_A(:,3), 360); %remove 2pi jumps and make continuous
euler_mah_A(:,3) = euler_mah_A(:,3) - euler_mah_A(StartIdx,3); %Mahony

warning('on','interpolation:interpolation:noextrap');

%% plotting
markerSize = 8; %pl

figPos = [105 105 645 660];

plotPos = [0.1 0.15 0.70 0.70];

plot21Pos = [0.1 0.53 0.86 0.40];
plot22Pos = plot21Pos + [0 -0.47 0 0];

plot1Pos = [0.1 0.69 0.86 0.25];
plot2Pos = plot1Pos + [0 -0.31 0 0];
plot3Pos = plot2Pos + [0 -0.31 0 0];

plot41Pos = [0.1 0.75 0.86 0.19];
plot42Pos = plot41Pos + [0 -0.23 0 0];
plot43Pos = plot42Pos + [0 -0.23 0 0];
plot44Pos = plot43Pos + [0 -0.23 0 0];

% Define locations from where tests start and stop (acceleration test and rotation test)

%ACC-test sequence: Linear motion with different accelerations
accStartIdx = find(time > 468, 1, 'first');
accStopIdx = find(time > 559, 1, 'first');
accStartIdx_IL = find(ILtime > 468, 1, 'first');
accStopIdx_IL = find(ILtime > 559, 1, 'first');
accStartIdx_SF = find(SFtime > 468, 1, 'first');
accStopIdx_SF = find(SFtime > 559, 1, 'first');
accStartIdx_K = find(Ktime > 468, 1, 'first');
accStopIdx_K = find(Ktime > 559, 1, 'first');

%sideways movement at x
accStartIdx_SF_x = find(SFtime > 468.5, 1, 'first');
accStopIdx_SF_x = find(SFtime > 496.5, 1, 'first');
%sideways movement at y
accStartIdx_SF_y = find(SFtime > 498, 1, 'first');
accStopIdx_SF_y = find(SFtime > 526, 1, 'first');
%up and down movement at z
accStartIdx_SF_z = find(SFtime > 530, 1, 'first');
accStopIdx_SF_z = find(SFtime > 558, 1, 'first');


%IMU-test sequence: free 6D motion
rotStartIdx = find(time > 564, 1, 'first');
rotStopIdx = find(time > 770, 1, 'first');

%raw accelerations
figIdx = 1;
figure(figIdx); clf; 

arrowOffset = 2;
arrowOffset_y = 0.15;

subplot(4,1,1); 
minY = -1;
maxY = 1;
plot(Ktime, vel_x, 'b', Ktime, vel_y, 'g--', Ktime, vel_z, 'r:', 'LineWidth', 1); 
set(gca, 'Position', plot41Pos, 'FontSize', 12, 'FontName', 'Times');
title('Reference measurement and accelerations during the acceleration test', 'FontSize', 14, 'FontName', 'Times');
legend('x', 'y', 'z', 'Location', 'NorthEast');
ylabel('Reference velocity (m/s)', 'FontSize', 12, 'FontName', 'Times');
axis([Ktime(accStartIdx_K) Ktime(accStopIdx_K) minY maxY]);
hold on;
plot(SFtime([accStartIdx_SF_x accStartIdx_SF_x]), [minY maxY], 'k:', SFtime([accStopIdx_SF_x accStopIdx_SF_x]), [minY maxY], 'k:', 'LineWidth', 1);
plot(SFtime([accStartIdx_SF_y accStartIdx_SF_y]), [minY maxY], 'k:', SFtime([accStopIdx_SF_y accStopIdx_SF_y]), [minY maxY], 'k:', 'LineWidth', 1);
plot(SFtime([accStartIdx_SF_z accStartIdx_SF_z]), [minY maxY], 'k:', SFtime([accStopIdx_SF_z accStopIdx_SF_z]), [minY maxY], 'k:', 'LineWidth', 1);


subplot(4,1,2); 
plot(ILtime(accStartIdx_IL), acc2_x(accStartIdx_IL), 'b--', ...
    SFtime(accStartIdx_SF:accStopIdx_SF), acc1_x(accStartIdx_SF:accStopIdx_SF), 'c', ...
    ILtime(accStartIdx_IL:accStopIdx_IL), acc2_x(accStartIdx_IL:accStopIdx_IL), 'b--'); 
set(gca, 'Position', plot42Pos, 'FontSize', 12, 'FontName', 'Times');
%title('Accelerations during the acceleration test', 'FontSize', 14, 'FontName', 'Times');
legend('A (Inertia-Link)', 'B (SparkFun)', 'Location', legendLocation(acc2_x(accStartIdx_IL:accStopIdx_IL)));
ylabel('acc_x (m/s^2)', 'FontSize', 12, 'FontName', 'Times');
hold on;

maxY = max([max(acc1_x(accStartIdx_SF:accStopIdx_SF)), max(acc2_x(accStartIdx_IL:accStopIdx_IL))]);
%maxY = max([1.1*maxY 0.9*maxY]);
minY = min([min(acc1_x(accStartIdx_SF:accStopIdx_SF)), min(acc2_x(accStartIdx_IL:accStopIdx_IL))]);
%minY = min([1.1*minY 0.9*minY]);
plot(SFtime([accStartIdx_SF_x accStartIdx_SF_x]), [minY maxY], 'k:', SFtime([accStopIdx_SF_x accStopIdx_SF_x]), [minY maxY], 'k:', 'LineWidth', 1);
plot(SFtime([accStartIdx_SF_y accStartIdx_SF_y]), [minY maxY], 'k:', SFtime([accStopIdx_SF_y accStopIdx_SF_y]), [minY maxY], 'k:', 'LineWidth', 1);
plot(SFtime([accStartIdx_SF_z accStartIdx_SF_z]), [minY maxY], 'k:', SFtime([accStopIdx_SF_z accStopIdx_SF_z]), [minY maxY], 'k:', 'LineWidth', 1);

axis([SFtime(accStartIdx_SF) SFtime(accStopIdx_SF) minY maxY]);

subplot(4,1,3); 
plot(SFtime(accStartIdx_SF:accStopIdx_SF), acc1_y(accStartIdx_SF:accStopIdx_SF), 'c', ...
    ILtime(accStartIdx_IL:accStopIdx_IL), acc2_y(accStartIdx_IL:accStopIdx_IL), 'b--'); 
set(gca, 'Position', plot43Pos, 'FontSize', 12, 'FontName', 'Times');
ylabel('acc_y (m/s^2)', 'FontSize', 12, 'FontName', 'Times');
hold on;

maxY = max([max(acc1_y(accStartIdx_SF:accStopIdx_SF)), max(acc2_y(accStartIdx_IL:accStopIdx_IL))]);
%maxY = max([1.1*maxY 0.9*maxY]);
minY = min([min(acc1_y(accStartIdx_SF:accStopIdx_SF)), min(acc2_y(accStartIdx_IL:accStopIdx_IL))]);
%minY = min([1.1*minY 0.9*minY]);
plot(SFtime([accStartIdx_SF_x accStartIdx_SF_x]), [minY maxY], 'k:', SFtime([accStopIdx_SF_x accStopIdx_SF_x]), [minY maxY], 'k:', 'LineWidth', 1);
plot(SFtime([accStartIdx_SF_y accStartIdx_SF_y]), [minY maxY], 'k:', SFtime([accStopIdx_SF_y accStopIdx_SF_y]), [minY maxY], 'k:', 'LineWidth', 1);
plot(SFtime([accStartIdx_SF_z accStartIdx_SF_z]), [minY maxY], 'k:', SFtime([accStopIdx_SF_z accStopIdx_SF_z]), [minY maxY], 'k:', 'LineWidth', 1);
axis([SFtime(accStartIdx_SF) SFtime(accStopIdx_SF) minY maxY]);

subplot(4,1,4); 
plot(SFtime(accStartIdx_SF:accStopIdx_SF), acc1_z(accStartIdx_SF:accStopIdx_SF), 'c', ...
    ILtime(accStartIdx_IL:accStopIdx_IL), acc2_z(accStartIdx_IL:accStopIdx_IL), 'b--'); 
set(gca, 'Position', plot44Pos, 'FontSize', 12, 'FontName', 'Times');
ylabel('acc_z (m/s^2)', 'FontSize', 12, 'FontName', 'Times');
xlabel('time (s)', 'FontSize', 12, 'FontName', 'Times');
hold on;

maxY = max([max(acc1_z(accStartIdx_SF:accStopIdx_SF)), max(acc2_z(accStartIdx_IL:accStopIdx_IL))]);
%maxY = max([1.1*maxY 0.9*maxY]);
minY = min([min(acc1_z(accStartIdx_SF:accStopIdx_SF)), min(acc2_z(accStartIdx_IL:accStopIdx_IL))]);
%minY = min([1.1*minY 0.9*minY]);
plot(SFtime([accStartIdx_SF_x accStartIdx_SF_x]), [minY maxY], 'k:', SFtime([accStopIdx_SF_x accStopIdx_SF_x]), [minY maxY], 'k:', 'LineWidth', 1);
plot(SFtime([accStartIdx_SF_y accStartIdx_SF_y]), [minY maxY], 'k:', SFtime([accStopIdx_SF_y accStopIdx_SF_y]), [minY maxY], 'k:', 'LineWidth', 1);
plot(SFtime([accStartIdx_SF_z accStartIdx_SF_z]), [minY maxY], 'k:', SFtime([accStopIdx_SF_z accStopIdx_SF_z]), [minY maxY], 'k:', 'LineWidth', 1);

plot([SFtime(accStartIdx_SF_x)+arrowOffset SFtime(accStopIdx_SF_x)-arrowOffset], [minY minY]+arrowOffset_y, 'k:', ...
    SFtime(accStartIdx_SF_x)+arrowOffset, minY+arrowOffset_y, 'k<', SFtime(accStopIdx_SF_x)-arrowOffset, minY+arrowOffset_y, 'k>', 'LineWidth', 1, 'MarkerSize', markerSize);
text(SFtime(accStartIdx_SF_x)+2*arrowOffset, minY-arrowOffset_y, 'x-axis movement', 'FontSize', 12, 'FontName', 'Times');

plot([SFtime(accStartIdx_SF_y)+arrowOffset SFtime(accStopIdx_SF_y)-arrowOffset], [minY minY]+arrowOffset_y, 'k:', ...
    SFtime(accStartIdx_SF_y)+arrowOffset, minY+arrowOffset_y, 'k<', SFtime(accStopIdx_SF_y)-arrowOffset, minY+arrowOffset_y, 'k>', 'LineWidth', 1, 'MarkerSize', markerSize);
text(SFtime(accStartIdx_SF_y)+2*arrowOffset, minY-arrowOffset_y, 'y-axis movement', 'FontSize', 12, 'FontName', 'Times');

plot([SFtime(accStartIdx_SF_z)+arrowOffset SFtime(accStopIdx_SF_z)-arrowOffset], [minY minY]+arrowOffset_y, 'k:', ...
    SFtime(accStartIdx_SF_z)+arrowOffset, minY+arrowOffset_y, 'k<', SFtime(accStopIdx_SF_z)-arrowOffset, minY+arrowOffset_y, 'k>', 'LineWidth', 1, 'MarkerSize', markerSize);
text(SFtime(accStartIdx_SF_z)+2*arrowOffset, minY-arrowOffset_y, 'z-axis movement', 'FontSize', 12, 'FontName', 'Times');

axis([SFtime(accStartIdx_SF) SFtime(accStopIdx_SF) minY-4*arrowOffset_y maxY]);

set(gcf, 'Position', figPos + [(figIdx-1)*100 0 0 150]);
set(gcf, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperPosition', [0.63 0.63 19.72 28.41]);
set(gcf,'PaperPositionMode','auto');
saveas(gcf, ['figures/test_fig_' int2str(figIdx) '.fig']);
print('-depsc','-tiff','-r300',['figures/test_fig_' int2str(figIdx) '.eps'])
figIdx = figIdx + 1;


% DCM states
figure(figIdx); clf; 
subplot(3,1,1); 
plot(ILtime, x_hist2(:,1), 'b--',  SFtime, x_hist1(:,1), 'c', ILtime, -data.InertiaLink.M13, 'm-.', Ktime, R31, 'k--'); 
legend('DCM_A', 'DCM_B', 'Inertia-Link', 'Kuka', 'Location', legendLocation([x_hist2(:,1) -data.InertiaLink.M13]));
set(gca, 'Position', plot1Pos, 'FontSize', 12, 'FontName', 'Times');
xlabel('time (s)', 'FontSize', 12, 'FontName', 'Times');
ylabel('x1 (R31)', 'FontSize', 12, 'FontName', 'Times');
title('DCM states', 'FontSize', 14, 'FontName', 'Times');
axis([0 commonDuration -1.1 1.1]);

subplot(3,1,2); 
plot(ILtime, x_hist2(:,2), 'b--',  SFtime, x_hist1(:,2), 'c', ILtime, -data.InertiaLink.M23, 'm-.', Ktime, R32, 'k--'); 
set(gca, 'Position', plot2Pos, 'FontSize', 12, 'FontName', 'Times');
ylabel('x2 (R32)', 'FontSize', 12, 'FontName', 'Times');
axis([0 commonDuration -1.1 1.1]);

subplot(3,1,3); 
plot(ILtime, x_hist2(:,3), 'b--', SFtime, x_hist1(:,3), 'c', ILtime, -data.InertiaLink.M33, 'm-.', Ktime, R33, 'k--'); 
set(gca, 'Position', plot3Pos, 'FontSize', 12, 'FontName', 'Times');
xlabel('time (s)', 'FontSize', 12, 'FontName', 'Times');
ylabel('x3 (R33)', 'FontSize', 12, 'FontName', 'Times');
axis([0 commonDuration -1.1 1.1]);

set(gcf, 'Position', figPos + [(figIdx-1)*100 0 0 0]);
set(gcf, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperPosition', [0.63 0.63 19.72 28.41]);
set(gcf,'PaperPositionMode','auto');
saveas(gcf, ['figures/test_fig_' int2str(figIdx) '.fig']);
print('-depsc','-tiff','-r300',['figures/test_fig_' int2str(figIdx) '.eps'])
figIdx = figIdx + 1;


% yaw
disp(' ');
disp('Yaw errors:')
dt = 1/150;
yaw_DCM_error = (ypr_hist2(StartIdx:StopIdx,1)-yaw_in_IL)*180/pi;
yaw_DCM_delta_error = [0; yaw_DCM_error(2:end) - yaw_DCM_error(1:end-1)] ./dt;
yaw_DCM_error_acc = yaw_DCM_error(accStartIdx:accStopIdx)-yaw_DCM_error(accStartIdx);
yaw_DCM_error_rot = yaw_DCM_error(rotStartIdx:rotStopIdx)-yaw_DCM_error(rotStartIdx);
yaw_DCM_RMSE_acc = sqrt(mean(yaw_DCM_error_acc.^2))
yaw_DCM_RMSE_rot = sqrt(mean(yaw_DCM_error_rot.^2))

yaw_sf_error = yaw_sf_in_IL - yaw_in_IL*180/pi;
yaw_sf_delta_error = [0; yaw_sf_error(2:end) - yaw_sf_error(1:end-1)] ./dt;
yaw_sf_error_acc = yaw_sf_error(accStartIdx:accStopIdx)-yaw_sf_error(accStartIdx);
yaw_sf_error_rot = yaw_sf_error(rotStartIdx:rotStopIdx)-yaw_sf_error(rotStartIdx);
yaw_sf_RMSE_acc = sqrt(mean(yaw_sf_error_acc.^2))
yaw_sf_RMSE_rot = sqrt(mean(yaw_sf_error_rot.^2))

yaw_mad_B_error = euler_mad_B(StartIdx:StopIdx,3) - yaw_in_IL*180/pi;
yaw_mad_B_delta_error = [0; yaw_mad_B_error(2:end) - yaw_mad_B_error(1:end-1)] ./dt;
yaw_mad_B_error_acc = yaw_mad_B_error(accStartIdx:accStopIdx)-yaw_mad_B_error(accStartIdx);
yaw_mad_B_error_rot = yaw_mad_B_error(rotStartIdx:rotStopIdx)-yaw_mad_B_error(rotStartIdx);
yaw_mad_B_RMSE_acc = sqrt(mean(yaw_mad_B_error_acc.^2))
yaw_mad_B_RMSE_rot = sqrt(mean(yaw_mad_B_error_rot.^2))

yaw_mah_B_error = euler_mah_B(StartIdx:StopIdx,3) - yaw_in_IL*180/pi;
yaw_mah_B_delta_error = [0; yaw_mah_B_error(2:end) - yaw_mah_B_error(1:end-1)] ./dt;
yaw_mah_B_error_acc = yaw_mah_B_error(accStartIdx:accStopIdx)-yaw_mah_B_error(accStartIdx);
yaw_mah_B_error_rot = yaw_mah_B_error(rotStartIdx:rotStopIdx)-yaw_mah_B_error(rotStartIdx);
yaw_mah_B_RMSE_acc = sqrt(mean(yaw_mah_B_error_acc.^2))
yaw_mah_B_RMSE_rot = sqrt(mean(yaw_mah_B_error_rot.^2))

yaw_mad_A_error = euler_mad_A(StartIdx:StopIdx,3) - yaw_in_IL*180/pi;
yaw_mad_A_delta_error = [0; yaw_mad_A_error(2:end) - yaw_mad_A_error(1:end-1)] ./dt;
yaw_mad_A_error_acc = yaw_mad_A_error(accStartIdx:accStopIdx)-yaw_mad_A_error(accStartIdx);
yaw_mad_A_error_rot = yaw_mad_A_error(rotStartIdx:rotStopIdx)-yaw_mad_A_error(rotStartIdx);
yaw_mad_A_RMSE_acc = sqrt(mean(yaw_mad_A_error_acc.^2))
yaw_mad_A_RMSE_rot = sqrt(mean(yaw_mad_A_error_rot.^2))

yaw_mah_A_error = euler_mah_A(StartIdx:StopIdx,3) - yaw_in_IL*180/pi;
yaw_mah_A_delta_error = [0; yaw_mah_A_error(2:end) - yaw_mah_A_error(1:end-1)] ./dt;
yaw_mah_A_error_acc = yaw_mah_A_error(accStartIdx:accStopIdx)-yaw_mah_A_error(accStartIdx);
yaw_mah_A_error_rot = yaw_mah_A_error(rotStartIdx:rotStopIdx)-yaw_mah_A_error(rotStartIdx);
yaw_mah_A_RMSE_acc = sqrt(mean(yaw_mah_A_error_acc.^2))
yaw_mah_A_RMSE_rot = sqrt(mean(yaw_mah_A_error_rot.^2))

yaw_IL_error = (yaw_IL(StartIdx:StopIdx)-yaw_in_IL)*180/pi;
yaw_IL_delta_error = [0; yaw_IL_error(2:end) - yaw_IL_error(1:end-1)] ./dt;
yaw_IL_error_acc = yaw_IL_error(accStartIdx:accStopIdx)-yaw_IL_error(accStartIdx);
yaw_IL_error_rot = yaw_IL_error(rotStartIdx:rotStopIdx)-yaw_IL_error(rotStartIdx);
yaw_IL_RMSE_acc = sqrt(mean(yaw_IL_error_acc.^2))
yaw_IL_RMSE_rot = sqrt(mean(yaw_IL_error_rot.^2))

%pitch
disp(' ');
disp('Pitch errors:')
pitch_DCM_error = (ypr_hist2(StartIdx:StopIdx,2)-pitch_in_IL)*180/pi;
pitch_DCM_RMSE_acc = sqrt(mean(pitch_DCM_error(accStartIdx:accStopIdx).^2))
pitch_DCM_RMSE_rot = sqrt(mean(pitch_DCM_error(rotStartIdx:rotStopIdx).^2))

pitch_sf_error = pitch_sf_in_IL - pitch_in_IL*180/pi;
pitch_sf_RMSE_acc = sqrt(mean(pitch_sf_error(accStartIdx:accStopIdx).^2))
pitch_sf_RMSE_rot = sqrt(mean(pitch_sf_error(rotStartIdx:rotStopIdx).^2))

pitch_mad_B_error = euler_mad_B(StartIdx:StopIdx,2) - pitch_in_IL*180/pi;
pitch_mad_B_RMSE_acc = sqrt(mean(pitch_mad_B_error(accStartIdx:accStopIdx).^2))
pitch_mad_B_RMSE_rot = sqrt(mean(pitch_mad_B_error(rotStartIdx:rotStopIdx).^2))

pitch_mah_B_error = euler_mah_B(StartIdx:StopIdx,2) - pitch_in_IL*180/pi;
pitch_mah_B_RMSE_acc = sqrt(mean(pitch_mah_B_error(accStartIdx:accStopIdx).^2))
pitch_mah_B_RMSE_rot = sqrt(mean(pitch_mah_B_error(rotStartIdx:rotStopIdx).^2))

pitch_mad_A_error = euler_mad_A(StartIdx:StopIdx,2) - pitch_in_IL*180/pi;
pitch_mad_A_RMSE_acc = sqrt(mean(pitch_mad_A_error(accStartIdx:accStopIdx).^2))
pitch_mad_A_RMSE_rot = sqrt(mean(pitch_mad_A_error(rotStartIdx:rotStopIdx).^2))

pitch_mah_A_error = euler_mah_A(StartIdx:StopIdx,2) - pitch_in_IL*180/pi;
pitch_mah_A_RMSE_acc = sqrt(mean(pitch_mah_A_error(accStartIdx:accStopIdx).^2))
pitch_mah_A_RMSE_rot = sqrt(mean(pitch_mah_A_error(rotStartIdx:rotStopIdx).^2))

pitch_IL_error = (pitch_IL(StartIdx:StopIdx)-pitch_in_IL)*180/pi;
pitch_IL_RMSE_acc = sqrt(mean(pitch_IL_error(accStartIdx:accStopIdx).^2))
pitch_IL_RMSE_rot = sqrt(mean(pitch_IL_error(rotStartIdx:rotStopIdx).^2))

%roll
disp(' ');
disp('Roll errors:')
roll_DCM_error = (ypr_hist2(StartIdx:StopIdx,3)-roll_in_IL)*180/pi;
roll_DCM_RMSE_acc = sqrt(mean(roll_DCM_error(accStartIdx:accStopIdx).^2))
roll_DCM_RMSE_rot = sqrt(mean(roll_DCM_error(rotStartIdx:rotStopIdx).^2))

roll_sf_error = roll_sf_in_IL - roll_in_IL*180/pi;
roll_sf_RMSE_acc = sqrt(mean(roll_sf_error(accStartIdx:accStopIdx).^2))
roll_sf_RMSE_rot = sqrt(mean(roll_sf_error(rotStartIdx:rotStopIdx).^2))

roll_mad_B_error = euler_mad_B(StartIdx:StopIdx,1) - roll_in_IL*180/pi;
roll_mad_B_RMSE_acc = sqrt(mean(roll_mad_B_error(accStartIdx:accStopIdx).^2))
roll_mad_B_RMSE_rot = sqrt(mean(roll_mad_B_error(rotStartIdx:rotStopIdx).^2))

roll_mah_B_error = euler_mah_B(StartIdx:StopIdx,1) - roll_in_IL*180/pi;
roll_mah_B_RMSE_acc = sqrt(mean(roll_mah_B_error(accStartIdx:accStopIdx).^2))
roll_mah_B_RMSE_rot = sqrt(mean(roll_mah_B_error(rotStartIdx:rotStopIdx).^2))

roll_mad_A_error = euler_mad_A(StartIdx:StopIdx,1) - roll_in_IL*180/pi;
roll_mad_A_RMSE_acc = sqrt(mean(roll_mad_A_error(accStartIdx:accStopIdx).^2))
roll_mad_A_RMSE_rot = sqrt(mean(roll_mad_A_error(rotStartIdx:rotStopIdx).^2))

roll_mah_A_error = euler_mah_A(StartIdx:StopIdx,1) - roll_in_IL*180/pi;
roll_mah_A_RMSE_acc = sqrt(mean(roll_mah_A_error(accStartIdx:accStopIdx).^2))
roll_mah_A_RMSE_rot = sqrt(mean(roll_mah_A_error(rotStartIdx:rotStopIdx).^2))

roll_IL_error = (roll_IL(StartIdx:StopIdx)-roll_in_IL)*180/pi;
roll_IL_RMSE_acc = sqrt(mean(roll_IL_error(accStartIdx:accStopIdx).^2))
roll_IL_RMSE_rot = sqrt(mean(roll_IL_error(rotStartIdx:rotStopIdx).^2))

disp(' ');
disp('Error matrices:')
RMSE_acc = [yaw_DCM_RMSE_acc yaw_sf_RMSE_acc yaw_mad_A_RMSE_acc yaw_mad_B_RMSE_acc yaw_mah_A_RMSE_acc yaw_mah_B_RMSE_acc yaw_IL_RMSE_acc; ...
    pitch_DCM_RMSE_acc pitch_sf_RMSE_acc pitch_mad_A_RMSE_acc pitch_mad_B_RMSE_acc pitch_mah_A_RMSE_acc pitch_mah_B_RMSE_acc pitch_IL_RMSE_acc; ...
    roll_DCM_RMSE_acc roll_sf_RMSE_acc roll_mad_A_RMSE_acc roll_mad_B_RMSE_acc roll_mah_A_RMSE_acc roll_mah_B_RMSE_acc roll_IL_RMSE_acc];

disp('RMSEs of the acceleration test');
textCell = [{'DCM_A', 'DCM_B', 'Madgwick_A', 'Madgwick_B', 'Mahony_A', 'Mahony_B', 'Inertia-Link'}; cell(size(RMSE_acc))];
for i = 1:size(RMSE_acc,1)
    for j = 1:size(RMSE_acc,2)
        textCell{i+1,j} = num2str(RMSE_acc(i,j), '%.2f');
    end
end
disp(textCell');
    
RMSE_rot = [yaw_DCM_RMSE_rot yaw_sf_RMSE_rot yaw_mad_A_RMSE_rot yaw_mad_B_RMSE_rot yaw_mah_A_RMSE_rot yaw_mah_B_RMSE_rot yaw_IL_RMSE_rot; ...
    pitch_DCM_RMSE_rot pitch_sf_RMSE_rot pitch_mad_A_RMSE_rot pitch_mad_B_RMSE_rot pitch_mah_A_RMSE_rot pitch_mah_B_RMSE_rot pitch_IL_RMSE_rot; ...
    roll_DCM_RMSE_rot roll_sf_RMSE_rot roll_mad_A_RMSE_rot roll_mad_B_RMSE_rot roll_mah_A_RMSE_rot roll_mah_B_RMSE_rot roll_IL_RMSE_rot];

disp('RMSEs of the rotation test');
textCell = [{'DCM_A', 'DCM_B', 'Madgwick_A', 'Madgwick_B', 'Mahony_A', 'Mahony_B', 'Inertia-Link'}; cell(size(RMSE_rot))];
for i = 1:size(RMSE_rot,1)
    for j = 1:size(RMSE_rot,2)
        textCell{i+1,j} = num2str(RMSE_rot(i,j), '%.2f');
    end
end
disp(textCell');


% box and whiskers plot for the acceleration test
groups_yaw = {'DCM_A', 'Madgwick_A', 'Mahony_A', 'DCM_B','Madgwick_B', 'Mahony_B', 'Inertia-Link_ '};

groups_rollpitch = cell(1, 2*size(groups_yaw,2));
groups_rollpitch(1:2:end) = groups_yaw;
groups_rollpitch(2:2:end) = groups_yaw;

figure(figIdx); clf; 
subplot(2,1,1);
boxplot([yaw_DCM_error_acc, yaw_mad_A_error_acc, yaw_mah_A_error_acc, yaw_sf_error_acc, yaw_mad_B_error_acc, yaw_mah_B_error_acc, yaw_IL_error_acc], ...
    groups_yaw, 'whisker', 15);
txt = findobj(gca,'Type','text');
set(txt, 'FontSize', 12, 'FontName', 'Times', 'Interpreter', 'tex', 'VerticalAlignment', 'middle');
set(gca, 'Position', plot21Pos, 'FontSize', 12, 'FontName', 'Times');
title('Error statistics in the acceleration test', 'FontSize', 14, 'FontName', 'Times');
ylabel('yaw errors (deg)');

subplot(2,1,2);

boxplot([pitch_DCM_error(accStartIdx:accStopIdx), roll_DCM_error(accStartIdx:accStopIdx), ...
    pitch_mad_A_error(accStartIdx:accStopIdx), roll_mad_A_error(accStartIdx:accStopIdx), ...
    pitch_mah_A_error(accStartIdx:accStopIdx), roll_mah_A_error(accStartIdx:accStopIdx), ...
    pitch_sf_error(accStartIdx:accStopIdx), roll_sf_error(accStartIdx:accStopIdx), ...
    pitch_mad_B_error(accStartIdx:accStopIdx), roll_mad_B_error(accStartIdx:accStopIdx), ...
    pitch_mah_B_error(accStartIdx:accStopIdx), roll_mah_B_error(accStartIdx:accStopIdx), ...
    pitch_IL_error(accStartIdx:accStopIdx), roll_IL_error(accStartIdx:accStopIdx)], ...
    groups_rollpitch, 'whisker', 15);

txt = findobj(gca,'Type','text');
set(txt, 'FontSize', 12, 'FontName', 'Times', 'Interpreter', 'tex', 'VerticalAlignment', 'middle');

set(gca, 'Position', plot22Pos, 'FontSize', 12, 'FontName', 'Times');
ylabel('combined roll and pitch errors (deg)');

set(gcf, 'Position', figPos + [(figIdx-1)*100 0 0 0]);
set(gcf, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperPosition', [0.63 0.63 19.72 28.41]);
set(gcf,'PaperPositionMode','auto');
saveas(gcf, ['figures/test_fig_' int2str(figIdx) '.fig']);
print('-depsc','-tiff','-r300',['figures/test_fig_' int2str(figIdx) '.eps'])
figIdx = figIdx + 1;   


% box and whiskers plot for the rotation test
figure(figIdx); clf; 
subplot(2,1,1);
boxplot([yaw_DCM_error_rot, yaw_mad_A_error_rot, yaw_mah_A_error_rot, yaw_sf_error_rot, yaw_mad_B_error_rot, yaw_mah_B_error_rot, yaw_IL_error_rot], ...
    groups_yaw, 'whisker', 15);
txt = findobj(gca,'Type','text');
set(txt, 'FontSize', 12, 'FontName', 'Times', 'Interpreter', 'tex', 'VerticalAlignment', 'middle');

set(gca, 'Position', plot21Pos, 'FontSize', 12, 'FontName', 'Times');
title('Error statistics in the rotation test', 'FontSize', 14, 'FontName', 'Times');
ylabel('yaw errors (deg)');

subplot(2,1,2);
boxplot([pitch_DCM_error(rotStartIdx:rotStopIdx), roll_DCM_error(rotStartIdx:rotStopIdx), ...
    pitch_mad_A_error(rotStartIdx:rotStopIdx), roll_mad_A_error(rotStartIdx:rotStopIdx), ...
    pitch_mah_A_error(rotStartIdx:rotStopIdx), roll_mah_A_error(rotStartIdx:rotStopIdx), ...
    pitch_sf_error(rotStartIdx:rotStopIdx), roll_sf_error(rotStartIdx:rotStopIdx), ...
    pitch_mad_B_error(rotStartIdx:rotStopIdx), roll_mad_B_error(rotStartIdx:rotStopIdx), ...
    pitch_mah_B_error(rotStartIdx:rotStopIdx), roll_mah_B_error(rotStartIdx:rotStopIdx), ...
    pitch_IL_error(rotStartIdx:rotStopIdx), roll_IL_error(rotStartIdx:rotStopIdx)], ...
    groups_rollpitch, 'whisker', 15);
txt = findobj(gca,'Type','text');
set(txt, 'FontSize', 12, 'FontName', 'Times', 'Interpreter', 'tex', 'VerticalAlignment', 'middle');

set(gca, 'Position', plot22Pos, 'FontSize', 12, 'FontName', 'Times');
ylabel('combined roll and pitch errors (deg)');

set(gcf, 'Position', figPos + [(figIdx-1)*100 0 0 0]);
set(gcf, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperPosition', [0.63 0.63 19.72 28.41]);
set(gcf,'PaperPositionMode','auto');
saveas(gcf, ['figures/test_fig_' int2str(figIdx) '.fig']);
print('-depsc','-tiff','-r300',['figures/test_fig_' int2str(figIdx) '.eps'])
figIdx = figIdx + 1;    


% correlation tests for delays
yaw_ref = yaw_in_IL*180/pi;
pitch_ref = pitch_in_IL*180/pi;
roll_ref = roll_in_IL*180/pi;

yaw_dcm = ypr_hist2(StartIdx:StopIdx,1)*180/pi;
pitch_dcm = ypr_hist2(StartIdx:StopIdx,2)*180/pi;
roll_dcm = ypr_hist2(StartIdx:StopIdx,3)*180/pi;

madgwick = euler_mad_A(StartIdx:StopIdx,:);
mahony = euler_mah_A(StartIdx:StopIdx,:);

yaw_il = yaw_IL(StartIdx:StopIdx)*180/pi;
pitch_il = pitch_IL(StartIdx:StopIdx)*180/pi;
roll_il = roll_IL(StartIdx:StopIdx)*180/pi;

n_lags = 30;
[c_dcm_yaw,lags] = rootMeanSquaredErrors(yaw_dcm, yaw_ref, n_lags);
[min_c_dcm_yaw, min_c_dcm_yaw_idx] = min(c_dcm_yaw);
c_madgwick_yaw = rootMeanSquaredErrors(madgwick(:,3), yaw_ref, n_lags);
[min_c_madgwick_yaw, min_c_madgwick_yaw_idx] = min(c_madgwick_yaw);
c_mahony_yaw = rootMeanSquaredErrors(mahony(:,3), yaw_ref, n_lags);
[min_c_mahony_yaw, min_c_mahony_yaw_idx] = min(c_mahony_yaw);
c_il_yaw = rootMeanSquaredErrors(yaw_il, yaw_ref, n_lags);
[min_c_il_yaw, min_c_il_yaw_idx] = min(c_il_yaw);

c_dcm_pitch = rootMeanSquaredErrors(pitch_dcm, pitch_ref, n_lags);
[min_c_dcm_pitch, min_c_dcm_pitch_idx] = min(c_dcm_pitch);
c_madgwick_pitch = rootMeanSquaredErrors(madgwick(:,2), pitch_ref, n_lags);
[min_c_madgwick_pitch, min_c_madgwick_pitch_idx] = min(c_madgwick_pitch);
c_mahony_pitch = rootMeanSquaredErrors(mahony(:,2), pitch_ref, n_lags);
[min_c_mahony_pitch, min_c_mahony_pitch_idx] = min(c_mahony_pitch);
c_il_pitch = rootMeanSquaredErrors(pitch_il, pitch_ref, n_lags);
[min_c_il_pitch, min_c_il_pitch_idx] = min(c_il_pitch);

c_dcm_roll = rootMeanSquaredErrors(roll_dcm, roll_ref, n_lags);
[min_c_dcm_roll, min_c_dcm_roll_idx] = min(c_dcm_roll);
c_madgwick_roll = rootMeanSquaredErrors(madgwick(:,1), roll_ref, n_lags);
[min_c_madgwick_roll, min_c_madgwick_roll_idx] = min(c_madgwick_roll);
c_mahony_roll = rootMeanSquaredErrors(mahony(:,1), roll_ref, n_lags);
[min_c_mahony_roll, min_c_mahony_roll_idx] = min(c_mahony_roll);
c_il_roll = rootMeanSquaredErrors(roll_il, roll_ref, n_lags);
[min_c_il_roll, min_c_il_roll_idx] = min(c_il_roll);

lags = lags*dt;


% mean squared errors plot
figure(figIdx); clf; 
subplot(3,1,1); 
plot(lags, c_dcm_yaw, 'b', lags, c_madgwick_yaw, 'g--', lags, c_mahony_yaw, 'r:', lags, c_il_yaw, 'm-.', 'LineWidth', 1, 'MarkerSize', markerSize);
set(gca, 'Position', plot1Pos, 'FontSize', 12, 'FontName', 'Times');
title('RMSEs to the reference measurement as a function of time delay', 'FontSize', 14, 'FontName', 'Times');
ylabel('yaw RMSE (deg)');
legend('DCM', 'Madgwick', 'Mahony', 'Inertia-Link', 'Location', legendLocation([c_dcm_yaw' c_madgwick_yaw' c_mahony_yaw' c_il_yaw']));
hold on;
plot(lags(min_c_dcm_yaw_idx), c_dcm_yaw(min_c_dcm_yaw_idx), 'b*', ...
    lags(min_c_madgwick_yaw_idx), c_madgwick_yaw(min_c_madgwick_yaw_idx), 'g^', ...
    lags(min_c_mahony_yaw_idx), c_mahony_yaw(min_c_mahony_yaw_idx), 'ro', ...
    lags(min_c_il_yaw_idx), c_il_yaw(min_c_il_yaw_idx), 'mp', 'LineWidth', 1, 'MarkerSize', markerSize);


subplot(3,1,2); 
plot(lags, c_dcm_pitch, 'b', lags, c_madgwick_pitch, 'g--', lags, c_mahony_pitch, 'r:', lags, c_il_pitch, 'm-.', 'LineWidth', 1, 'MarkerSize', markerSize);
set(gca, 'Position', plot2Pos, 'FontSize', 12, 'FontName', 'Times');
ylabel('pitch RMSE (deg)');
hold on;
plot(lags(min_c_dcm_pitch_idx), c_dcm_pitch(min_c_dcm_pitch_idx), 'b*', ...
    lags(min_c_madgwick_pitch_idx), c_madgwick_pitch(min_c_madgwick_pitch_idx), 'g^', ...
    lags(min_c_mahony_pitch_idx), c_mahony_pitch(min_c_mahony_pitch_idx), 'ro', ...
    lags(min_c_il_pitch_idx), c_il_pitch(min_c_il_pitch_idx), 'mp', 'LineWidth', 1, 'MarkerSize', markerSize);


subplot(3,1,3); 
plot(lags, c_dcm_roll, 'b', lags, c_madgwick_roll, 'g--', lags, c_mahony_roll, 'r:', lags, c_il_roll, 'm-.', 'LineWidth', 1, 'MarkerSize', markerSize);
set(gca, 'Position', plot3Pos, 'FontSize', 12, 'FontName', 'Times');
ylabel('roll RMSE (deg)');
xlabel('time delay to reference (s)');
hold on;
plot(lags(min_c_dcm_roll_idx), c_dcm_roll(min_c_dcm_roll_idx), 'b*', ...
    lags(min_c_madgwick_roll_idx), c_madgwick_roll(min_c_madgwick_roll_idx), 'g^', ...
    lags(min_c_mahony_roll_idx), c_mahony_roll(min_c_mahony_roll_idx), 'ro', ...
    lags(min_c_il_roll_idx), c_il_roll(min_c_il_roll_idx), 'mp', 'LineWidth', 1, 'MarkerSize', markerSize);


set(gcf, 'Position', figPos + [(figIdx-1)*100 0 0 0]);
set(gcf, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperPosition', [0.63 0.63 19.72 28.41]);
set(gcf,'PaperPositionMode','auto');
saveas(gcf, ['figures/test_fig_' int2str(figIdx) '.fig']);
print('-depsc','-tiff','-r300',['figures/test_fig_' int2str(figIdx) '.eps'])
figIdx = figIdx + 1;    


% Yaw, pitch and roll plot
figure(figIdx); clf; 
subplot(3,1,1); 
plot(time, yaw_dcm, 'b', time, madgwick(:,3), 'g--', time, mahony(:,3), 'r:', time, yaw_il, 'm-.', time, yaw_ref, 'k--', 'LineWidth', 1); 
set(gca, 'Position', plot1Pos, 'FontSize', 12, 'FontName', 'Times');
legend('DCM', 'Madgwick', 'Mahony', 'Inertia-Link', 'Reference', 'Location', legendLocation([yaw_dcm yaw_ref madgwick(:,3) mahony(:,3)]));
hold on;
maxY = max([max(yaw_dcm), max(madgwick(:,3)), max(mahony(:,3)), max(yaw_il), max(yaw_ref)]);
maxY = max([1.1*maxY 0.9*maxY]);
minY = min([min(yaw_dcm), min(madgwick(:,3)), min(mahony(:,3)), min(yaw_il), min(yaw_ref)]);
minY = min([1.1*minY 0.9*minY]);
plot(time([accStartIdx accStartIdx]), [minY maxY], 'k:', time([accStopIdx accStopIdx]), [minY maxY], 'k:', 'LineWidth', 1);
plot(time([rotStartIdx rotStartIdx]), [minY maxY], 'k:', time([rotStopIdx rotStopIdx]), [minY maxY], 'k:', 'LineWidth', 1);
ylabel('yaw (deg)', 'FontSize', 12, 'FontName', 'Times');
title('Euler angles', 'FontSize', 14, 'FontName', 'Times');
axis([time(1) time(end) minY maxY]);

subplot(3,1,2); 
plot(time, pitch_dcm, 'b', time, madgwick(:,2), 'g--', time, mahony(:,2), 'r:', time, pitch_il, 'm-.', time, pitch_ref, 'k--', 'LineWidth', 1);
hold on;
maxY = max([max(pitch_dcm), max(madgwick(:,2)), max(mahony(:,2)), max(pitch_il), max(pitch_ref)]);
maxY = max([1.1*maxY 0.9*maxY]);
minY = min([min(pitch_dcm), min(madgwick(:,2)), min(mahony(:,2)), min(pitch_il), min(pitch_ref)]);
minY = min([1.1*minY 0.9*minY]);
plot(time([accStartIdx accStartIdx]), [minY maxY], 'k:', time([accStopIdx accStopIdx]), [minY maxY], 'k:', 'LineWidth', 1);
plot(time([rotStartIdx rotStartIdx]), [minY maxY], 'k:', time([rotStopIdx rotStopIdx]), [minY maxY], 'k:', 'LineWidth', 1);
set(gca, 'Position', plot2Pos, 'FontSize', 12, 'FontName', 'Times');
ylabel('pitch (deg)', 'FontSize', 12, 'FontName', 'Times');
axis([time(1) time(end) minY maxY]);

subplot(3,1,3); 
plot(time, roll_dcm, 'b', time, madgwick(:,1), 'g--', time, mahony(:,1), 'r:', time, roll_il, 'm-.', time, roll_ref, 'k--', 'LineWidth', 1); 
set(gca, 'Position', plot3Pos, 'FontSize', 12, 'FontName', 'Times');
hold on;
maxY = max([max(roll_dcm), max(madgwick(:,1)), max(mahony(:,1)), max(roll_il), max(roll_ref)]);
maxY = max([1.1*maxY 0.9*maxY]);
minY = min([min(roll_dcm), min(madgwick(:,1)), min(mahony(:,1)), min(roll_il), min(roll_ref)]);
minY = min([1.1*minY 0.9*minY]);
arrowOffset = 5;
plot(time([accStartIdx accStartIdx]), [minY maxY], 'k:', time([accStopIdx accStopIdx]), [minY maxY], 'k:', 'LineWidth', 1);
plot([time(accStartIdx)+arrowOffset time(accStopIdx)-arrowOffset], [minY minY]+arrowOffset, 'k:', ...
    time(accStartIdx)+arrowOffset, minY+arrowOffset, 'k<', time(accStopIdx)-arrowOffset, minY+arrowOffset, 'k>', 'LineWidth', 1, 'MarkerSize', markerSize);
text(time(accStartIdx)+2*arrowOffset, minY-arrowOffset, 'Acceleration test', 'FontSize', 12, 'FontName', 'Times');
plot(time([rotStartIdx rotStartIdx]), [minY maxY], 'k:', time([rotStopIdx rotStopIdx]), [minY maxY], 'k:', 'LineWidth', 1);
plot([time(rotStartIdx)+arrowOffset time(rotStopIdx)-arrowOffset], [minY minY]+arrowOffset, 'k:', ...
    time(rotStartIdx)+arrowOffset, minY+arrowOffset, 'k<', time(rotStopIdx)-arrowOffset, minY+arrowOffset, 'k>', 'LineWidth', 1, 'MarkerSize', markerSize);
text(time(rotStartIdx)+2*arrowOffset, minY-arrowOffset, 'Rotation test', 'FontSize', 12, 'FontName', 'Times');
xlabel('time (s)', 'FontSize', 12, 'FontName', 'Times');
ylabel('roll (deg)', 'FontSize', 12, 'FontName', 'Times');

axis([time(1) time(end) minY-20 maxY]);

set(gcf, 'Position', figPos + [(figIdx-1)*100 0 0 0]);
set(gcf, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperPosition', [0.63 0.63 19.72 28.41]);
set(gcf,'PaperPositionMode','auto');
saveas(gcf, ['figures/test_fig_' int2str(figIdx) '.fig']);
print('-depsc','-tiff','-r300',['figures/test_fig_' int2str(figIdx) '.eps'])
figIdx = figIdx + 1;    


% error plot with reference
figure(figIdx); clf; 
subplot(4,1,1);
plot(time, yaw_ref, 'b', time, pitch_ref, 'g--', time, roll_ref, 'r:', 'LineWidth', 1); 
set(gca, 'Position', plot41Pos, 'FontSize', 12, 'FontName', 'Times');
legend('yaw', 'pitch', 'roll', 'Location', legendLocation([yaw_ref pitch_ref roll_ref]));
title('The reference measurement and measurement errors', 'FontSize', 14, 'FontName', 'Times');
ylabel('reference angles (deg)', 'FontSize', 12, 'FontName', 'Times');

hold on;
maxY = max([max(yaw_ref), max(pitch_ref), max(roll_ref)]);
maxY = max([1.1*maxY 0.9*maxY]);
minY = min([min(yaw_ref), min(pitch_ref), min(roll_ref)]);
minY = min([1.1*minY 0.9*minY]);
arrowOffset = 5;
plot(time([accStartIdx accStartIdx]), [minY maxY], 'k:', time([accStopIdx accStopIdx]), [minY maxY], 'k:', 'LineWidth', 1);
plot(time([rotStartIdx rotStartIdx]), [minY maxY], 'k:', time([rotStopIdx rotStopIdx]), [minY maxY], 'k:', 'LineWidth', 1);

axis([time(1) time(end) minY maxY]);


subplot(4,1,2); 
plot(time, yaw_DCM_error, 'b', time, yaw_mad_A_error, 'g--', time, yaw_mah_A_error, 'r:', time, yaw_IL_error, 'm-.', 'LineWidth', 1, 'EraseMode', 'xor'); 
set(gca, 'Position', plot42Pos, 'FontSize', 12, 'FontName', 'Times');
legend('DCM', 'Madgwick', 'Mahony', 'Inertia-Link', 'Location', legendLocation([yaw_DCM_error yaw_IL_error yaw_mad_A_error yaw_mah_A_error]));
hold on;
ylabel('yaw errors (deg)', 'FontSize', 12, 'FontName', 'Times');
set(gca, 'XLim', [time(1) time(end)]);

maxY = max([max(yaw_DCM_error), max(yaw_mad_A_error), max(yaw_mah_A_error), max(yaw_IL_error)]);
maxY = max([1.1*maxY 0.9*maxY]);
minY = min([min(yaw_DCM_error), min(yaw_mad_A_error), min(yaw_mah_A_error), min(yaw_IL_error)]);
minY = min([1.1*minY 0.9*minY]);

plot(time([accStartIdx accStartIdx]), [minY maxY], 'k:', time([accStopIdx accStopIdx]), [minY maxY], 'k:', 'LineWidth', 1);
plot(time([rotStartIdx rotStartIdx]), [minY maxY], 'k:', time([rotStopIdx rotStopIdx]), [minY maxY], 'k:', 'LineWidth', 1);

axis([time(1) time(end) minY maxY]);


subplot(4,1,3); 
plot(time, pitch_DCM_error, 'b', time, pitch_mad_A_error, 'g--', time, pitch_mah_A_error, 'r:', time, pitch_IL_error, 'm-.', 'LineWidth', 1, 'EraseMode', 'xor'); 
hold on;
set(gca, 'Position', plot43Pos, 'FontSize', 12, 'FontName', 'Times');
ylabel('pitch errors (deg)', 'FontSize', 12, 'FontName', 'Times');

maxY = max([max(pitch_DCM_error), max(pitch_mad_A_error), max(pitch_mah_A_error), max(pitch_IL_error)]);
maxY = max([1.1*maxY 0.9*maxY]);
minY = min([min(pitch_DCM_error), min(pitch_mad_A_error), min(pitch_mah_A_error), min(pitch_IL_error)]);
minY = min([1.1*minY 0.9*minY]);

plot(time([accStartIdx accStartIdx]), [minY maxY], 'k:', time([accStopIdx accStopIdx]), [minY maxY], 'k:', 'LineWidth', 1);
plot(time([rotStartIdx rotStartIdx]), [minY maxY], 'k:', time([rotStopIdx rotStopIdx]), [minY maxY], 'k:', 'LineWidth', 1);

axis([time(1) time(end) minY maxY]);


subplot(4,1,4); 
plot(time, roll_DCM_error, 'b', time, roll_mad_A_error, 'g--', time, roll_mah_A_error, 'r:', time, roll_IL_error, 'm-.', 'LineWidth', 1, 'EraseMode', 'xor');
hold on;
set(gca, 'Position', plot44Pos, 'FontSize', 12, 'FontName', 'Times');
xlabel('time (s)', 'FontSize', 12, 'FontName', 'Times');
ylabel('roll errors (deg)', 'FontSize', 12, 'FontName', 'Times');

maxY = max([max(roll_DCM_error), max(roll_mad_A_error), max(roll_mah_A_error), max(roll_IL_error)]);
maxY = max([1.1*maxY 0.9*maxY]);
minY = min([min(roll_DCM_error), min(roll_mad_A_error), min(roll_mah_A_error), min(roll_IL_error)]);
minY = min([1.1*minY 0.9*minY]);

arrowOffset_y = 0.5;

plot(time([accStartIdx accStartIdx]), [minY maxY], 'k:', time([accStopIdx accStopIdx]), [minY maxY], 'k:', 'LineWidth', 1);
plot([time(accStartIdx)+arrowOffset time(accStopIdx)-arrowOffset], [minY minY]+arrowOffset_y, 'k:', ...
    time(accStartIdx)+arrowOffset, minY+arrowOffset_y, 'k<', time(accStopIdx)-arrowOffset, minY+arrowOffset_y, 'k>', 'LineWidth', 1, 'MarkerSize', markerSize);
text(time(accStartIdx)+2*arrowOffset, minY-arrowOffset_y, 'Acceleration test', 'FontSize', 12, 'FontName', 'Times');
plot(time([rotStartIdx rotStartIdx]), [minY maxY], 'k:', time([rotStopIdx rotStopIdx]), [minY maxY], 'k:', 'LineWidth', 1);
plot([time(rotStartIdx)+arrowOffset time(rotStopIdx)-arrowOffset], [minY minY]+arrowOffset_y, 'k:', ...
    time(rotStartIdx)+arrowOffset, minY+arrowOffset_y, 'k<', time(rotStopIdx)-arrowOffset, minY+arrowOffset_y, 'k>', 'LineWidth', 1, 'MarkerSize', markerSize);
text(time(rotStartIdx)+2*arrowOffset, minY-arrowOffset_y, 'Rotation test', 'FontSize', 12, 'FontName', 'Times');

axis([time(1) time(end) minY-4*arrowOffset_y maxY]);

set(gcf, 'Position', figPos + [(figIdx-1)*100 0 0 150]);
set(gcf, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperPosition', [0.63 0.63 19.72 28.41]);
set(gcf,'PaperPositionMode','auto');
saveas(gcf, ['figures/test_fig_' int2str(figIdx) '.fig']);
print('-depsc','-tiff','-r300',['figures/test_fig_' int2str(figIdx) '.eps'])
figIdx = figIdx + 1;


% bias plot (states with variances)
figure(figIdx); clf; 

x_sigma = sqrt(P_diag_hist2(StartIdx:StopIdx,4));
y_sigma = sqrt(P_diag_hist2(StartIdx:StopIdx,5));
z_sigma = sqrt(P_diag_hist2(StartIdx:StopIdx,6));

plotLims = 1.25*[mean(x_sigma), mean(y_sigma), mean(z_sigma)];

bias_x = x_hist2(StartIdx:StopIdx,4)*180/pi;
bias_y = x_hist2(StartIdx:StopIdx,5)*180/pi;
bias_z = x_hist2(StartIdx:StopIdx,6)*180/pi;

subplot(3,1,1);
plot([time(1) time(end)], [addedGyroBias addedGyroBias]*180/pi, 'k--', ...
    time, bias_x, 'b', time, (addedGyroBias-x_sigma)*180/pi, 'b-.', ...
    time, (addedGyroBias+x_sigma)*180/pi, 'b-.', 'LineWidth', 1, 'MarkerSize', markerSize); 
h = legend('Reference', 'Estimate', '1-sigma');
set(gca, 'Position', plot1Pos, 'FontSize', 12, 'FontName', 'Times');

ylabel('x_b_i_a_s (deg/s)', 'FontSize', 12, 'FontName', 'Times');
title('Bias estimates and 1-sigma distances of standard deviations', 'FontSize', 14, 'FontName', 'Times');
set(gca, 'XLim', [time(1) time(end)], 'YLim', [addedGyroBias-plotLims(1) addedGyroBias+plotLims(1)]*180/pi);

subplot(3,1,2);
plot([time(1) time(end)], [addedGyroBias addedGyroBias]*180/pi, 'k--', ...
    time, bias_y, 'b', time, (addedGyroBias-y_sigma)*180/pi, 'b-.', ...
    time, (addedGyroBias+y_sigma)*180/pi, 'b-.', 'LineWidth', 1, 'MarkerSize', markerSize); 
set(gca, 'Position', plot2Pos, 'FontSize', 12, 'FontName', 'Times');

ylabel('y_b_i_a_s (deg/s)', 'FontSize', 12, 'FontName', 'Times');
set(gca, 'XLim', [time(1) time(end)], 'YLim', [addedGyroBias-plotLims(2) addedGyroBias+plotLims(2)]*180/pi);

subplot(3,1,3);
plot([time(1) time(end)], [addedGyroBias addedGyroBias]*180/pi, 'k--', ...
    time, bias_z, 'b', time, (addedGyroBias-z_sigma)*180/pi, 'b-.', ...
    time, (addedGyroBias+z_sigma)*180/pi, 'b-.', 'LineWidth', 1, 'MarkerSize', markerSize); 
set(gca, 'Position', plot3Pos, 'FontSize', 12, 'FontName', 'Times');

xlabel('time (s)', 'FontSize', 12, 'FontName', 'Times');
ylabel('z_b_i_a_s (deg/s)', 'FontSize', 12, 'FontName', 'Times');
set(gca, 'XLim', [time(1) time(end)], 'YLim', [addedGyroBias-plotLims(3) addedGyroBias+plotLims(3)]*180/pi);

set(gcf, 'Position', figPos + [(figIdx-1)*100 0 0 0]);
set(gcf, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperPosition', [0.63 0.63 19.72 28.41]);
set(gcf,'PaperPositionMode','auto');
saveas(gcf, ['figures/test_fig_' int2str(figIdx) '.fig']);
print('-depsc','-tiff','-r300',['figures/test_fig_' int2str(figIdx) '.eps'])




