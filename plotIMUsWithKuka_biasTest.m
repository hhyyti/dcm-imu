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

%% Configuration!
% select to test biases using Acceleration test (true) or Rotation test (false)
ACCtest = true; 

% use c/compile.m before setting this to true (tested only with linux)
use_C_versions = false; 

% fetch GPL licensed algortihms in matlab and C from here before use: 
% http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/ 
% Add the C code to c folder and matlab functions in this folder 
% (please add also the quaternion library used by Madgwick).
% I apologize this mess. The Madgwick's GPS licenced code is separated from 
% this MIT licenced product to avoid infecting this code with GPL licence.
use_comparisonAlgorithms = false; 


if (ACCtest)
    disp('Calculating and plotting IMU bias test for Acceleration data');
else
    disp('Calculating and plotting IMU bias test for Rotation data');
end
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


if (ACCtest)
    %ACC-test sequence: Linear motion with different accelerations
    StartIdx = find(ILtime > 468, 1, 'first');
    StopIdx = find(ILtime > 559, 1, 'first');
else
    %IMU-test sequence: free 6D motion
    StartIdx = find(ILtime > 564, 1, 'first');
    StopIdx = find(ILtime > 770, 1, 'first');
end

testName = 'Rotation test';
if (ACCtest)
    testName = 'Acceleration test';
end

markerSize = 8;

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

%% getting reference trajectory
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

%% calibrate IMU data
imusWithKukaCalibration;


if (~use_comparisonAlgorithms)
    disp('To get reference algorithms, download their implementations from here:');
    disp('http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/');
    disp(' ');
end

%% start bias test computation (compute everything with test region using different biases)
% test and computation is done only for InertiaLink data
testBiases = ((0:0.5:7)*pi/180)';
%testBiases = ((0:0.5:1)*pi/180)';

yaw_DCM_RMSE = zeros(size(testBiases));
yaw_e1_RMSE = zeros(size(testBiases));
yaw_e2_RMSE = zeros(size(testBiases));

pitch_DCM_RMSE = zeros(size(testBiases));
pitch_e1_RMSE = zeros(size(testBiases));
pitch_e2_RMSE = zeros(size(testBiases));

roll_DCM_RMSE = zeros(size(testBiases));
roll_e1_RMSE = zeros(size(testBiases));
roll_e2_RMSE = zeros(size(testBiases));

acc2 = [acc2_x(StartIdx:StopIdx), acc2_y(StartIdx:StopIdx), acc2_z(StartIdx:StopIdx)];
gyro2_noBias = [w2_x(StartIdx:StopIdx), w2_y(StartIdx:StopIdx), w2_z(StartIdx:StopIdx)];
time = ILtime(StartIdx:StopIdx);
deltaT = [0; time(2:end) - time(1:end-1)];

%resample Kuka reference data to InertiaLink time
warning('off','interpolation:interpolation:noextrap');
yaw_in_IL = resample(timeseries(yaw, Ktime),time);
yaw_in_IL = yaw_in_IL.Data;
yaw_in_IL = yaw_in_IL - yaw_in_IL(1); %starting yaw from zero (same as IMU estimates)
pitch_in_IL = resample(timeseries(pitch, Ktime),time);
pitch_in_IL = pitch_in_IL.Data;
roll_in_IL = resample(timeseries(roll, Ktime),time);
roll_in_IL = roll_in_IL.Data;
warning('on','interpolation:interpolation:noextrap');

disp('Starting to compute RMSEs using different biases');

%check if matlabpool is opened
scheduler = findResource();
if (matlabpool('size') < scheduler.ClusterSize)
    if (matlabpool('size') > 0)
        matlabpool close;
    end
        
    matlabpool(scheduler.ClusterSize);
end

warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
parfor i = 2:(length(testBiases))
    disp(['iteration ' int2str(i) ', bias ' num2str(testBiases(i)*180/pi) ' deg/s']);
    
    % add constant bias for testing purposes
    addedGyroBias = testBiases(i);
    gyro2 = gyro2_noBias + addedGyroBias*ones(size(gyro2_noBias));
  
    if (use_C_versions) 
        [x_hist2, ypr_hist2, a_hist2, P_diag_hist2] = ...
            DCM_IMU_C(gyro2, acc2, deltaT);
        
        if (use_comparisonAlgorithms)
            quaternion1 = Madgwick_IMU_C(gyro2, acc2);
            quaternion2 = Mahony_IMU_C(gyro2, acc2);        
        else
            quaternion1 = nan(size(gyro2,1),4);
            quaternion2 = nan(size(gyro2,1),4);
        end
    else
        DCM = DCM_IMU();
        
        if (use_comparisonAlgorithms)
            addpath('quaternion_library');      % include quaternion library
            
            IMU1t = MadgwickAHRS('SamplePeriod', 1/150, 'Beta', 0.1);
            IMU2t = MahonyAHRS('SamplePeriod', 1/150, 'Kp', 0.5);
        
            quaternion1 = zeros(length(time), 4);
            quaternion2 = zeros(length(time), 4);
        else
            quaternion1 = nan(length(time), 4);
            quaternion2 = nan(length(time), 4);
        end
        
        x_hist2 = zeros(length(time), 6);
        ypr_hist2 = zeros(length(time), 3);
        P_diag_hist2 = zeros(length(time), 6);
        for t = 1:length(time)
            DCM.UpdateIMU(gyro2(t,:), acc2(t,:), deltaT(t));	% gyroscope units must be radians
            x_hist2(t, :) = DCM.state';
            ypr_hist2(t, :) = [DCM.yaw, DCM.pitch, DCM.roll];
            P_diag_hist2(t, :) = diag(DCM.P)';

            if (use_comparisonAlgorithms)
                IMU1t.UpdateIMU(gyro2(t,:), acc2(t,:));	% gyroscope units must be radians
                quaternion1(t, :) = IMU1t.Quaternion;

                IMU2t.UpdateIMU(gyro2(t,:), acc2(t,:));	% gyroscope units must be radians
                quaternion2(t, :) = IMU2t.Quaternion;
            end
        end        
    end
  

    % Plot algorithm output as Euler angles
    % The first and third Euler angles in the sequence (phi and psi) become
    % unreliable when the middle angles of the sequence (theta) approaches ~90
    % degrees. This problem commonly referred to as Gimbal Lock.
    % See: http://en.wikipedia.org/wiki/Gimbal_lock

    % use conjugate for sensor frame relative to Earth and convert to degrees.
    if (use_comparisonAlgorithms)
        euler1 = quatern2euler(quaternConj(quaternion1)) * (180/pi);	
        euler2 = quatern2euler(quaternConj(quaternion2)) * (180/pi);	
    else
        euler1 = nan(size(quaternion1,1),3);
        euler2 = nan(size(quaternion2,1),3);
    end	

    % compute errors
    
    %yaw pitch roll errors (only InertiaLink data)
    ypr_hist2(:,1) = un_modulo(ypr_hist2(:,1), 2*pi); %remove 2pi jumps and make continuous
    ypr_hist2(:,1) = ypr_hist2(:,1) - ypr_hist2(1,1); %DCM

    euler1(:,3) = un_modulo(euler1(:,3), 360); %remove 2pi jumps and make continuous
    euler1(:,3) = euler1(:,3) - euler1(1,3); %Madgwick

    euler2(:,3) = un_modulo(euler2(:,3), 360); %remove 2pi jumps and make continuous
    euler2(:,3) = euler2(:,3) - euler2(1,3); %Mahony

    %yaw
    yaw_DCM_error = (ypr_hist2(:,1)-yaw_in_IL)*180/pi;
    yaw_DCM_RMSE(i) = sqrt(mean(yaw_DCM_error.^2));
    yaw_e1_error = euler1(:,3) - yaw_in_IL*180/pi;
    yaw_e1_RMSE(i) = sqrt(mean(yaw_e1_error.^2));
    yaw_e2_error = euler2(:,3) - yaw_in_IL*180/pi;
    yaw_e2_RMSE(i) = sqrt(mean(yaw_e2_error.^2));

    %pitch
    pitch_DCM_error = (ypr_hist2(:,2)-pitch_in_IL)*180/pi;
    pitch_DCM_RMSE(i) = sqrt(mean(pitch_DCM_error.^2));
    pitch_e1_error = euler1(:,2) - pitch_in_IL*180/pi;
    pitch_e1_RMSE(i) = sqrt(mean(pitch_e1_error.^2));
    pitch_e2_error = euler2(:,2) - pitch_in_IL*180/pi;
    pitch_e2_RMSE(i) = sqrt(mean(pitch_e2_error.^2));

    %roll
    roll_DCM_error = (ypr_hist2(:,3)-roll_in_IL)*180/pi;
    roll_DCM_RMSE(i) = sqrt(mean(roll_DCM_error.^2));
    roll_e1_error = euler1(:,1) - roll_in_IL*180/pi;
    roll_e1_RMSE(i) = sqrt(mean(roll_e1_error.^2));
    roll_e2_error = euler2(:,1) - roll_in_IL*180/pi;
    roll_e2_RMSE(i) = sqrt(mean(roll_e2_error.^2));
end
warning('on','MATLAB:mir_warning_maybe_uninitialized_temporary');
%% plotting

% compute same for smallest and some other index (for separate plotting)
figIdx = 1;
j_idx = [1 3];
for j = 1:2
    i = j_idx(j);
    
    disp(['iteration ' int2str(i) ', bias ' num2str(testBiases(i)*180/pi) ' deg/s']);
    
    % add constant bias for testing purposes
    addedGyroBias = testBiases(i);
    gyro2 = gyro2_noBias + addedGyroBias*ones(size(gyro2_noBias));

    if (use_C_versions) 
        [x_hist2, ypr_hist2, a_hist2, P_diag_hist2] = ...
            DCM_IMU_C(gyro2, acc2, deltaT);
        
        if (use_comparisonAlgorithms)
            quaternion1 = Madgwick_IMU_C(gyro2, acc2);
            quaternion2 = Mahony_IMU_C(gyro2, acc2);        
        else
            quaternion1 = nan(size(gyro2,1),4);
            quaternion2 = nan(size(gyro2,1),4);
        end 
    else
        DCM = DCM_IMU();
        
        if (use_comparisonAlgorithms)
            addpath('quaternion_library');      % include quaternion library

            IMU1 = MadgwickAHRS('SamplePeriod', 1/150, 'Beta', 0.1);
            IMU2 = MahonyAHRS('SamplePeriod', 1/150, 'Kp', 0.5);
        
            quaternion1 = zeros(length(time), 4);
            quaternion2 = zeros(length(time), 4);
        else
            quaternion1 = nan(length(time), 4);
            quaternion2 = nan(length(time), 4);
        end
        
        x_hist2 = zeros(length(time), 6);
        ypr_hist2 = zeros(length(time), 3);
        P_diag_hist2 = zeros(length(time), 6);
        for t = 1:length(time)
            DCM.UpdateIMU(gyro2(t,:), acc2(t,:), deltaT(t));	% gyroscope units must be radians
            x_hist2(t, :) = DCM.state';
            ypr_hist2(t, :) = [DCM.yaw, DCM.pitch, DCM.roll];
            P_diag_hist2(t, :) = diag(DCM.P)';

            if (use_comparisonAlgorithms)
                IMU1.UpdateIMU(gyro2(t,:), acc2(t,:));	% gyroscope units must be radians
                quaternion1(t, :) = IMU1.Quaternion;

                IMU2.UpdateIMU(gyro2(t,:), acc2(t,:));	% gyroscope units must be radians
                quaternion2(t, :) = IMU2.Quaternion;
            end
        end        
    end

    % Plot algorithm output as Euler angles
    % The first and third Euler angles in the sequence (phi and psi) become
    % unreliable when the middle angles of the sequence (theta) approaches ~90
    % degrees. This problem commonly referred to as Gimbal Lock.
    % See: http://en.wikipedia.org/wiki/Gimbal_lock

    % use conjugate for sensor frame relative to Earth and convert to degrees.
    if (use_comparisonAlgorithms)
        euler1 = quatern2euler(quaternConj(quaternion1)) * (180/pi);	
        euler2 = quatern2euler(quaternConj(quaternion2)) * (180/pi);	
    else
        euler1 = nan(size(quaternion1,1),3);
        euler2 = nan(size(quaternion2,1),3);
    end

    % compute errors

    % yaw pitch roll errors (only InertiaLink data)
    ypr_hist2(:,1) = un_modulo(ypr_hist2(:,1), 2*pi); %remove 2pi jumps and make continuous
    ypr_hist2(:,1) = ypr_hist2(:,1) - ypr_hist2(1,1); %DCM

    euler1(:,3) = un_modulo(euler1(:,3), 360); %remove 2pi jumps and make continuous
    euler1(:,3) = euler1(:,3) - euler1(1,3); %Madgwick

    euler2(:,3) = un_modulo(euler2(:,3), 360); %remove 2pi jumps and make continuous
    euler2(:,3) = euler2(:,3) - euler2(1,3); %Mahony

    %yaw
    yaw_DCM_error = (ypr_hist2(:,1)-yaw_in_IL)*180/pi;
    yaw_DCM_RMSE(i) = sqrt(mean(yaw_DCM_error.^2));
    yaw_e1_error = euler1(:,3) - yaw_in_IL*180/pi;
    yaw_e1_RMSE(i) = sqrt(mean(yaw_e1_error.^2));
    yaw_e2_error = euler2(:,3) - yaw_in_IL*180/pi;
    yaw_e2_RMSE(i) = sqrt(mean(yaw_e2_error.^2));

    %pitch
    pitch_DCM_error = (ypr_hist2(:,2)-pitch_in_IL)*180/pi;
    pitch_DCM_RMSE(i) = sqrt(mean(pitch_DCM_error.^2));
    pitch_e1_error = euler1(:,2) - pitch_in_IL*180/pi;
    pitch_e1_RMSE(i) = sqrt(mean(pitch_e1_error.^2));
    pitch_e2_error = euler2(:,2) - pitch_in_IL*180/pi;
    pitch_e2_RMSE(i) = sqrt(mean(pitch_e2_error.^2));

    %roll
    roll_DCM_error = (ypr_hist2(:,3)-roll_in_IL)*180/pi;
    roll_DCM_RMSE(i) = sqrt(mean(roll_DCM_error.^2));
    roll_e1_error = euler1(:,1) - roll_in_IL*180/pi;
    roll_e1_RMSE(i) = sqrt(mean(roll_e1_error.^2));
    roll_e2_error = euler2(:,1) - roll_in_IL*180/pi;
    roll_e2_RMSE(i) = sqrt(mean(roll_e2_error.^2));

    % plotting
    % Yaw, pitch and roll plot
    figure(figIdx); clf; 
    subplot(3,1,1); 
    yaw_dcm = ypr_hist2(:,1)*180/pi;
    pitch_dcm = ypr_hist2(:,2)*180/pi;
    roll_dcm = ypr_hist2(:,3)*180/pi;
    
    plot(time, yaw_dcm, 'b',  time, euler1(:,3), 'g--', time, euler2(:,3), 'r:', time, yaw_in_IL*180/pi, 'k--', 'LineWidth', 1); 
    set(gca, 'Position', plot1Pos, 'FontSize', 12, 'FontName', 'Times');
    legend('DCM', 'Madgwick', 'Mahony', 'Reference', 'Location', legendLocation([yaw_dcm euler1(:,3) euler2(:,3)]));

    ylabel('yaw (deg)', 'FontSize', 12, 'FontName', 'Times');
    title(['Euler angles with a bias of ' num2str(addedGyroBias*180/pi, '%.0f') ' deg/s'], 'FontSize', 14, 'FontName', 'Times');
    set(gca, 'XLim', [time(1) time(end)]);

    subplot(3,1,2); 
    plot(time, pitch_dcm, 'b',  time, euler1(:,2), 'g--', time, euler2(:,2), 'r:', time, pitch_in_IL*180/pi, 'k--', 'LineWidth', 1);
    set(gca, 'Position', plot2Pos, 'FontSize', 12, 'FontName', 'Times');
    ylabel('pitch (deg)', 'FontSize', 12, 'FontName', 'Times');
    set(gca, 'XLim', [time(1) time(end)]);

    subplot(3,1,3); 
    plot(time, roll_dcm, 'b',  time, euler1(:,1), 'g--', time, euler2(:,1), 'r:', time, roll_in_IL*180/pi, 'k--', 'LineWidth', 1); 
    set(gca, 'Position', plot3Pos, 'FontSize', 12, 'FontName', 'Times');
    xlabel('time (s)', 'FontSize', 12, 'FontName', 'Times');
    ylabel('roll (deg)', 'FontSize', 12, 'FontName', 'Times');
    set(gca, 'XLim', [time(1) time(end)]);

    set(gcf, 'Position', figPos + [(figIdx-1)*100 0 0 0]);
    set(gcf, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperPosition', [0.63 0.63 19.72 28.41]);
    set(gcf,'PaperPositionMode','auto');
    saveas(gcf, ['figures/bias_fig_' testName(1:3) '_' int2str(figIdx) '.fig']);
    print('-depsc','-tiff','-r300',['figures/bias_fig_' testName(1:3) '_' int2str(figIdx) '.eps'])
    figIdx = figIdx + 1;    
    
    % error plot
    figure(figIdx); clf; 

    subplot(3,1,1); 
    plot(time, yaw_DCM_error, 'b', time, yaw_e1_error, 'g--', time, yaw_e2_error, 'r:', 'LineWidth', 1); 
    set(gca, 'Position', plot1Pos, 'FontSize', 12, 'FontName', 'Times');
    legend('DCM', 'Madgwick', 'Mahony', 'Location', legendLocation([yaw_DCM_error yaw_e1_error yaw_e2_error]));
    ylabel('yaw error (deg)', 'FontSize', 12, 'FontName', 'Times');
    title(['Errors to the reference measurement with a bias of ' num2str(addedGyroBias*180/pi, '%.0f') ' deg/s'], 'FontSize', 14, 'FontName', 'Times');
    set(gca, 'XLim', [time(1) time(end)]);

    subplot(3,1,2); 
    plot(time, pitch_DCM_error, 'b', time, pitch_e1_error, 'g--', time, pitch_e2_error, 'r:', 'LineWidth', 1); 
    set(gca, 'Position', plot2Pos, 'FontSize', 12, 'FontName', 'Times');
    ylabel('pitch error (deg)', 'FontSize', 12, 'FontName', 'Times');
    set(gca, 'XLim', [time(1) time(end)]);

    subplot(3,1,3); 
    plot(time, roll_DCM_error, 'b',  time, roll_e1_error, 'g--', time, roll_e2_error, 'r:', 'LineWidth', 1);
    set(gca, 'Position', plot3Pos, 'FontSize', 12, 'FontName', 'Times');
    xlabel('time (s)', 'FontSize', 12, 'FontName', 'Times');
    ylabel('roll error (deg)', 'FontSize', 12, 'FontName', 'Times');
    set(gca, 'XLim', [time(1) time(end)]);

    set(gcf, 'Position', figPos + [(figIdx-1)*100 0 0 0]);
    set(gcf, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperPosition', [0.63 0.63 19.72 28.41]);
    set(gcf,'PaperPositionMode','auto');
    saveas(gcf, ['figures/bias_fig_' testName(1:3) '_' int2str(figIdx) '.fig']);
    print('-depsc','-tiff','-r300',['figures/bias_fig_' testName(1:3) '_' int2str(figIdx) '.eps'])
    figIdx = figIdx + 1;
    
    % bias plot
    figure(figIdx); clf; 

    x_sigma = sqrt(P_diag_hist2(:,4));
    y_sigma = sqrt(P_diag_hist2(:,5));
    z_sigma = sqrt(P_diag_hist2(:,6));

    bias_x = x_hist2(:,4)*180/pi;
    bias_y = x_hist2(:,5)*180/pi;
    bias_z = x_hist2(:,6)*180/pi;
    
    plotLims = 1.25*[mean(x_sigma(1:floor(end/3))), mean(y_sigma(1:floor(end/3))), mean(z_sigma(1:floor(end/3)))];

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
    saveas(gcf, ['figures/bias_fig_' testName(1:3) '_' int2str(figIdx) '.fig']);
    print('-depsc','-tiff','-r300',['figures/bias_fig_' testName(1:3) '_' int2str(figIdx) '.eps'])
    figIdx = figIdx + 1;
    
    pause(0.2);
end

%% statistics plot
figure(figIdx); clf; 

subplot(3,1,1); 
semilogy(testBiases*180/pi, yaw_DCM_RMSE, 'b*-', testBiases*180/pi, yaw_e1_RMSE, 'g^-', testBiases*180/pi, yaw_e2_RMSE, 'ro-', 'LineWidth', 1, 'MarkerSize', markerSize); 
set(gca, 'Position', plot1Pos, 'FontSize', 12, 'FontName', 'Times');
legend('DCM', 'Madgwick', 'Mahony', 'Location', 'SouthEast');
ylabel('yaw RMSE (deg)', 'FontSize', 12, 'FontName', 'Times');
title(['RMSEs as a function of gyroscope bias (' testName ')'], 'FontSize', 14, 'FontName', 'Times');

subplot(3,1,2); 
semilogy(testBiases*180/pi, pitch_DCM_RMSE, 'b*-', testBiases*180/pi, pitch_e1_RMSE, 'g^-', testBiases*180/pi, pitch_e2_RMSE, 'ro-', 'LineWidth', 1, 'MarkerSize', markerSize); 
set(gca, 'Position', plot2Pos, 'FontSize', 12, 'FontName', 'Times');
ylabel('pitch RMSE (deg)', 'FontSize', 12, 'FontName', 'Times');

subplot(3,1,3); 
semilogy(testBiases*180/pi, roll_DCM_RMSE, 'b*-', testBiases*180/pi, roll_e1_RMSE, 'g^-', testBiases*180/pi, roll_e2_RMSE, 'ro-', 'LineWidth', 1, 'MarkerSize', markerSize); 
set(gca, 'Position', plot3Pos, 'FontSize', 12, 'FontName', 'Times');
xlabel('Constant gyroscope bias (deg/s)', 'FontSize', 12, 'FontName', 'Times');
ylabel('roll RMSE (deg)', 'FontSize', 12, 'FontName', 'Times');

set(gcf, 'Position', figPos + [(figIdx-1)*100 0 0 0]);
set(gcf, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperPosition', [0.63 0.63 19.72 28.41]);
set(gcf,'PaperPositionMode','auto');
saveas(gcf, ['figures/bias_fig_' testName(1:3) '_' int2str(figIdx) '.fig']);
print('-depsc','-tiff','-r300',['figures/bias_fig_' testName(1:3) '_' int2str(figIdx) '.eps'])


