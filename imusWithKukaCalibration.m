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

% acceleration data calibration 1 is for IMU1 (IL), 2 for IMU2 (IL) 
calibTimes = [9, 23; 26, 40; 42, 56; 59, 73; 76, 89; 93, 107]' + 4;
angleCalibTimes = [116, 225; 234, 345; 351, 464]';

calibTimes2 = calibTimes + 763;
angleCalibTimes2 = angleCalibTimes + 763;

minCalibTime = min(calibTimes(:));
maxCalibTime = max(angleCalibTimes2(:));

ILdelay = 0.03; %seconds, delay measured in InertiaLink data
SFdelay = 0; %seconds, delay measured in SparkFun data

% calibration point 1

S1_x = zeros(6,1);
S1_y = zeros(6,1);
S1_z = zeros(6,1);
w1_Bx = zeros(6,1);
w1_By = zeros(6,1);
w1_Bz = zeros(6,1);
T_1 = zeros(6,1);

S2_x = zeros(6,1);
S2_y = zeros(6,1);
S2_z = zeros(6,1);
w2_Bx = zeros(6,1);
w2_By = zeros(6,1);
w2_Bz = zeros(6,1);

for i = 1:6
    SFrange = SFtime > calibTimes(1,i) & SFtime < calibTimes(2,i);
    ILrange = ILtime > calibTimes(1,i) & ILtime < calibTimes(2,i);

    S1_x(i) = mean(data.SparkFun6DOF.acc_x(SFrange)); 
    S1_y(i) = mean(data.SparkFun6DOF.acc_y(SFrange));
    S1_z(i) = mean(data.SparkFun6DOF.acc_z(SFrange));
    w1_Bx(i) = mean(data.SparkFun6DOF.w_x(SFrange));
    w1_By(i) = mean(data.SparkFun6DOF.w_y(SFrange));
    w1_Bz(i) = mean(data.SparkFun6DOF.w_z(SFrange));
    T_1(i) = mean(data.SparkFun6DOF.temperature(SFrange));
    
    S2_x(i) = mean(data.InertiaLink.acc_x(ILrange));
    S2_y(i) = mean(data.InertiaLink.acc_y(ILrange));
    S2_z(i) = mean(data.InertiaLink.acc_z(ILrange));
    w2_Bx(i) = mean(data.InertiaLink.w_x(ILrange));
    w2_By(i) = mean(data.InertiaLink.w_y(ILrange));
    w2_Bz(i) = mean(data.InertiaLink.w_z(ILrange));
end
T1_acc = mean(T_1); %average temperature during acceleration calibration

% initialize biases and gains 
% (assume that bias is near zero and gain near one)
B1 = zeros(3,1);
B2 = zeros(3,1);

G1 = ones(3,1);
G2 = ones(3,1);

% run iterative calibration algorithm for both Temperatures
for i = 1:1000
    [B1_new, G1_new] = IMUcalibIteration(S1_x, S1_y, S1_z, B1, G1);
    [B2_new, G2_new] = IMUcalibIteration(S2_x, S2_y, S2_z, B2, G2);
       
    Differences = abs(B1_new - B1) + abs(G1_new - G1) + ...
        abs(B2_new - B2) + abs(G2_new - G2);
    if (Differences < 0.000001) 
        break;
    end
    
    B1 = B1_new;
    G1 = G1_new;
    B2 = B2_new;
    G2 = G2_new;
end


% calibrate gyroscopes by comparing to KUKA data

% resample Kuka angular velocities to SFtime
minTime = min(angleCalibTimes(:));
maxTime = max(angleCalibTimes(:));
SFrange = SFtime > minTime & SFtime < maxTime;
T1_gyro = mean(data.SparkFun6DOF.temperature(SFrange));
angleCalibTime_SF = SFtime(SFrange);

w_x_in_SF = resample(timeseries(gyro(:,1), Ktime+SFdelay),angleCalibTime_SF);
w_x_in_SF = w_x_in_SF.Data;
w_y_in_SF = resample(timeseries(gyro(:,2), Ktime+SFdelay),angleCalibTime_SF);
w_y_in_SF = w_y_in_SF.Data;
w_z_in_SF = resample(timeseries(gyro(:,3), Ktime+SFdelay),angleCalibTime_SF);
w_z_in_SF = w_z_in_SF.Data;

w_x_error_SF = data.SparkFun6DOF.w_x(SFrange) - w_x_in_SF;
w_y_error_SF = data.SparkFun6DOF.w_y(SFrange) - w_y_in_SF;
w_z_error_SF = data.SparkFun6DOF.w_z(SFrange) - w_z_in_SF;

%filter out range where reference angular velocity is near zero
filt_dist = 0.0055;
x_filt_SF = (w_x_in_SF > filt_dist) | (w_x_in_SF < -filt_dist);
y_filt_SF = (w_y_in_SF > filt_dist) | (w_y_in_SF < -filt_dist);
z_filt_SF = (w_z_in_SF > filt_dist) | (w_z_in_SF < -filt_dist);

%resample Kuka angular velocities to ILtime
ILrange = ILtime > minTime & ILtime < maxTime;
angleCalibTime_IL = ILtime(ILrange);

w_x_in_IL = resample(timeseries(gyro(:,1), Ktime+ILdelay),angleCalibTime_IL);
w_x_in_IL = w_x_in_IL.Data;
w_y_in_IL = resample(timeseries(gyro(:,2), Ktime+ILdelay),angleCalibTime_IL);
w_y_in_IL = w_y_in_IL.Data;
w_z_in_IL = resample(timeseries(gyro(:,3), Ktime+ILdelay),angleCalibTime_IL);
w_z_in_IL = w_z_in_IL.Data;

w_x_error_IL = data.InertiaLink.w_x(ILrange) - w_x_in_IL;
w_y_error_IL = data.InertiaLink.w_y(ILrange) - w_y_in_IL;
w_z_error_IL = data.InertiaLink.w_z(ILrange) - w_z_in_IL;

%filter out range where reference angular velocity is near zero
x_filt_IL = (w_x_in_IL > filt_dist) | (w_x_in_IL < -filt_dist);
y_filt_IL = (w_y_in_IL > filt_dist) | (w_y_in_IL < -filt_dist);
z_filt_IL = (w_z_in_IL > filt_dist) | (w_z_in_IL < -filt_dist);

[w1_Gx, w1_Bx] = fitRegressionLine(w_x_in_SF(x_filt_SF), w_x_in_SF(x_filt_SF)+w_x_error_SF(x_filt_SF));
[w2_Gx, w2_Bx] = fitRegressionLine(w_x_in_IL(x_filt_IL), w_x_in_IL(x_filt_IL)+w_x_error_IL(x_filt_IL));
[w1_Gy, w1_By] = fitRegressionLine(w_y_in_SF(y_filt_SF), w_y_in_SF(y_filt_SF)+w_y_error_SF(y_filt_SF));
[w2_Gy, w2_By] = fitRegressionLine(w_y_in_IL(y_filt_IL), w_y_in_IL(y_filt_IL)+w_y_error_IL(y_filt_IL));
[w1_Gz, w1_Bz] = fitRegressionLine(w_z_in_SF(z_filt_SF), w_z_in_SF(z_filt_SF)+w_z_error_SF(z_filt_SF));
[w2_Gz, w2_Bz] = fitRegressionLine(w_z_in_IL(z_filt_IL), w_z_in_IL(z_filt_IL)+w_z_error_IL(z_filt_IL));


% calibration point 2

S1b_x = zeros(6,1);
S1b_y = zeros(6,1);
S1b_z = zeros(6,1);
w1b_Bx = zeros(6,1);
w1b_By = zeros(6,1);
w1b_Bz = zeros(6,1);
T_2 = zeros(6,1);

S2b_x = zeros(6,1);
S2b_y = zeros(6,1);
S2b_z = zeros(6,1);
w2b_Bx = zeros(6,1);
w2b_By = zeros(6,1);
w2b_Bz = zeros(6,1);

for i = 1:6
    SFrange = SFtime > calibTimes2(1,i) & SFtime < calibTimes2(2,i);
    ILrange = ILtime > calibTimes2(1,i) & ILtime < calibTimes2(2,i);

    S1b_x(i) = mean(data.SparkFun6DOF.acc_x(SFrange)); 
    S1b_y(i) = mean(data.SparkFun6DOF.acc_y(SFrange));
    S1b_z(i) = mean(data.SparkFun6DOF.acc_z(SFrange));
    w1b_Bx(i) = mean(data.SparkFun6DOF.w_x(SFrange));
    w1b_By(i) = mean(data.SparkFun6DOF.w_y(SFrange));
    w1b_Bz(i) = mean(data.SparkFun6DOF.w_z(SFrange));
    T_2(i) = mean(data.SparkFun6DOF.temperature(SFrange));
    
    S2b_x(i) = mean(data.InertiaLink.acc_x(ILrange));
    S2b_y(i) = mean(data.InertiaLink.acc_y(ILrange));
    S2b_z(i) = mean(data.InertiaLink.acc_z(ILrange));
    w2b_Bx(i) = mean(data.InertiaLink.w_x(ILrange));
    w2b_By(i) = mean(data.InertiaLink.w_y(ILrange));
    w2b_Bz(i) = mean(data.InertiaLink.w_z(ILrange));
end
T2_acc = mean(T_2); %average temperature during acceleration calibration

% initialize biases and gains 
% (assume that bias is near zero and gain near one)
B1b = zeros(3,1);
B2b = zeros(3,1);

G1b = ones(3,1);
G2b = ones(3,1);

% run iterative calibration algorithm for both Temperatures
for i = 1:1000
    [B1b_new, G1b_new] = IMUcalibIteration(S1b_x, S1b_y, S1b_z, B1b, G1b);
    [B2b_new, G2b_new] = IMUcalibIteration(S2b_x, S2b_y, S2b_z, B2b, G2b);
       
    Differences = abs(B1b_new - B1b) + abs(G1b_new - G1b) + ...
        abs(B2b_new - B2b) + abs(G2b_new - G2b);
    if (Differences < 0.000001) 
        break;
    end
    
    B1b = B1b_new;
    G1b = G1b_new;
    B2b = B2b_new;
    G2b = G2b_new;
end

% calibrate gyroscopes by comparing to KUKA data

% resample Kuka angular velocities to SFtime
minTime = min(angleCalibTimes2(:));
maxTime = max(angleCalibTimes2(:));
SFrange = SFtime > minTime & SFtime < maxTime;
T2_gyro = mean(data.SparkFun6DOF.temperature(SFrange));
angleCalibTime_SF = SFtime(SFrange);

w_x_in_SF = resample(timeseries(gyro(:,1), Ktime+SFdelay),angleCalibTime_SF);
w_x_in_SF = w_x_in_SF.Data;
w_y_in_SF = resample(timeseries(gyro(:,2), Ktime+SFdelay),angleCalibTime_SF);
w_y_in_SF = w_y_in_SF.Data;
w_z_in_SF = resample(timeseries(gyro(:,3), Ktime+SFdelay),angleCalibTime_SF);
w_z_in_SF = w_z_in_SF.Data;

w_x_error_SF = data.SparkFun6DOF.w_x(SFrange) - w_x_in_SF;
w_y_error_SF = data.SparkFun6DOF.w_y(SFrange) - w_y_in_SF;
w_z_error_SF = data.SparkFun6DOF.w_z(SFrange) - w_z_in_SF;

% filter out range where reference angular velocity is near zero
filt_dist = 0.0055;
x_filt_SF = (w_x_in_SF > filt_dist) | (w_x_in_SF < -filt_dist);
y_filt_SF = (w_y_in_SF > filt_dist) | (w_y_in_SF < -filt_dist);
z_filt_SF = (w_z_in_SF > filt_dist) | (w_z_in_SF < -filt_dist);

% resample Kuka angular velocities to ILtime
ILrange = ILtime > minTime & ILtime < maxTime;
angleCalibTime_IL = ILtime(ILrange);

w_x_in_IL = resample(timeseries(gyro(:,1), Ktime+ILdelay),angleCalibTime_IL);
w_x_in_IL = w_x_in_IL.Data;
w_y_in_IL = resample(timeseries(gyro(:,2), Ktime+ILdelay),angleCalibTime_IL);
w_y_in_IL = w_y_in_IL.Data;
w_z_in_IL = resample(timeseries(gyro(:,3), Ktime+ILdelay),angleCalibTime_IL);
w_z_in_IL = w_z_in_IL.Data;

w_x_error_IL = data.InertiaLink.w_x(ILrange) - w_x_in_IL;
w_y_error_IL = data.InertiaLink.w_y(ILrange) - w_y_in_IL;
w_z_error_IL = data.InertiaLink.w_z(ILrange) - w_z_in_IL;

% filter out range where reference angular velocity is near zero
x_filt_IL = (w_x_in_IL > filt_dist) | (w_x_in_IL < -filt_dist);
y_filt_IL = (w_y_in_IL > filt_dist) | (w_y_in_IL < -filt_dist);
z_filt_IL = (w_z_in_IL > filt_dist) | (w_z_in_IL < -filt_dist);

[w1b_Gx, w1b_Bx] = fitRegressionLine(w_x_in_SF(x_filt_SF), w_x_in_SF(x_filt_SF)+w_x_error_SF(x_filt_SF));
[w2b_Gx, w2b_Bx] = fitRegressionLine(w_x_in_IL(x_filt_IL), w_x_in_IL(x_filt_IL)+w_x_error_IL(x_filt_IL));
[w1b_Gy, w1b_By] = fitRegressionLine(w_y_in_SF(y_filt_SF), w_y_in_SF(y_filt_SF)+w_y_error_SF(y_filt_SF));
[w2b_Gy, w2b_By] = fitRegressionLine(w_y_in_IL(y_filt_IL), w_y_in_IL(y_filt_IL)+w_y_error_IL(y_filt_IL));
[w1b_Gz, w1b_Bz] = fitRegressionLine(w_z_in_SF(z_filt_SF), w_z_in_SF(z_filt_SF)+w_z_error_SF(z_filt_SF));
[w2b_Gz, w2b_Bz] = fitRegressionLine(w_z_in_IL(z_filt_IL), w_z_in_IL(z_filt_IL)+w_z_error_IL(z_filt_IL));

% accelerometer parameters for sparkfun
a1_b_Ax = (B1b(1) - B1(1)) / (T2_acc - T1_acc);
a1_b_Ay = (B1b(2) - B1(2)) / (T2_acc - T1_acc);
a1_b_Az = (B1b(3) - B1(3)) / (T2_acc - T1_acc);

a1_b_Bx = B1(1) - a1_b_Ax*T1_acc;
a1_b_By = B1(2) - a1_b_Ay*T1_acc;
a1_b_Bz = B1(3) - a1_b_Az*T1_acc;

a1_g_Ax = (G1b(1) - G1(1)) / (T2_acc - T1_acc);
a1_g_Ay = (G1b(2) - G1(2)) / (T2_acc - T1_acc);
a1_g_Az = (G1b(3) - G1(3)) / (T2_acc - T1_acc);

a1_g_Bx = G1(1) - a1_g_Ax*T1_acc;
a1_g_By = G1(2) - a1_g_Ay*T1_acc;
a1_g_Bz = G1(3) - a1_g_Az*T1_acc;

% gyroscope parameters for sparkfun
w1_b_Ax = (w1b_Bx - w1_Bx) / (T2_gyro - T1_gyro);
w1_b_Ay = (w1b_By - w1_By) / (T2_gyro - T1_gyro);
w1_b_Az = (w1b_Bz - w1_Bz) / (T2_gyro - T1_gyro);

w1_b_Bx = w1_Bx - w1_b_Ax*T1_gyro;
w1_b_By = w1_By - w1_b_Ay*T1_gyro;
w1_b_Bz = w1_Bz - w1_b_Az*T1_gyro;

w1_g_Ax = (w1b_Gx - w1_Gx) / (T2_gyro - T1_gyro);
w1_g_Ay = (w1b_Gy - w1_Gy) / (T2_gyro - T1_gyro);
w1_g_Az = (w1b_Gz - w1_Gz) / (T2_gyro - T1_gyro);

w1_g_Bx = w1_Gx - w1_g_Ax*T1_gyro;
w1_g_By = w1_Gy - w1_g_Ay*T1_gyro;
w1_g_Bz = w1_Gz - w1_g_Az*T1_gyro;


% print output for matlab
disp(' ');
disp('Simple calibration model for gyroscopes and accelerometers');
disp('model: calibrated = (measurement - bias)/gain');
disp(' ');
disp('Sparkfun 6DOF digital:');
disp(['w1_bias_xyz = [' num2str(w1_Bx, '%.8f') ...
    ', ' num2str(w1_By, '%.8f') ', ' num2str(w1_Bz, '%.8f') '];']);
disp(['w1_gain_xyz = [' num2str(w1_Gx, '%.8f') ...
    ', ' num2str(w1_Gy, '%.8f') ', ' num2str(w1_Gz, '%.8f') '];']);
disp(['a1_bias_xyz = [' num2str(B1(1), '%.8f') ...
    ', ' num2str(B1(2), '%.8f') ', ' num2str(B1(3), '%.8f') '];']);
disp(['a1_gain_xyz = [' num2str(G1(1), '%.8f') ...
    ', ' num2str(G1(2), '%.8f') ', ' num2str(G1(3), '%.8f') '];']);
disp(' ');
disp('InertiaLink IMU:');
disp(['w2_bias_xyz = [' num2str(w2_Bx, '%.8f') ...
    ', ' num2str(w2_By, '%.8f') ', ' num2str(w2_Bz, '%.8f') '];']);
disp(['w2_gain_xyz = [' num2str(w2_Gx, '%.8f') ...
    ', ' num2str(w2_Gy, '%.8f') ', ' num2str(w2_Gz, '%.8f') '];']);
disp(['a2_bias_xyz = [' num2str(B2(1), '%.8f') ...
    ', ' num2str(B2(2), '%.8f') ', ' num2str(B2(3), '%.8f') '];']);
disp(['a2_gain_xyz = [' num2str(G2(1), '%.8f') ...
    ', ' num2str(G2(2), '%.8f') ', ' num2str(G2(3), '%.8f') '];']);
disp(' ');


disp(' ');
disp('Temperature based calibration model for gyroscopes and accelerometers');
disp('model: calibrated = (measurement - bias(T))/gain(T)');
disp('          bias(T) = bias_A*T + bias_B');
disp('          gain(T) = gain_A*T + gain_B');
disp(' ');
disp('Sparkfun 6DOF digital (only SparkFun has temperature data):');
disp(['w1_bias_a_xyz = [' num2str(w1_b_Ax, '%.8f') ...
    ', ' num2str(w1_b_Ay, '%.8f') ', ' num2str(w1_b_Az, '%.8f') '];']);
disp(['w1_bias_b_xyz = [' num2str(w1_b_Bx, '%.8f') ...
    ', ' num2str(w1_b_By, '%.8f') ', ' num2str(w1_b_Bz, '%.8f') '];']);

disp(['w1_gain_a_xyz = [' num2str(w1_g_Ax, '%.8f') ...
    ', ' num2str(w1_g_Ay, '%.8f') ', ' num2str(w1_g_Az, '%.8f') '];']);
disp(['w1_gain_b_xyz = [' num2str(w1_g_Bx, '%.8f') ...
    ', ' num2str(w1_g_By, '%.8f') ', ' num2str(w1_g_Bz, '%.8f') '];']);

disp(['a1_bias_a_xyz = [' num2str(a1_b_Ax, '%.8f') ...
    ', ' num2str(a1_b_Ay, '%.8f') ', ' num2str(a1_b_Az, '%.8f') '];']);
disp(['a1_bias_b_xyz = [' num2str(a1_b_Bx, '%.8f') ...
    ', ' num2str(a1_b_By, '%.8f') ', ' num2str(a1_b_Bz, '%.8f') '];']);

disp(['a1_gain_a_xyz = [' num2str(a1_g_Ax, '%.8f') ...
    ', ' num2str(a1_g_Ay, '%.8f') ', ' num2str(a1_g_Az, '%.8f') '];']);
disp(['a1_gain_b_xyz = [' num2str(a1_g_Bx, '%.8f') ...
    ', ' num2str(a1_g_By, '%.8f') ', ' num2str(a1_g_Bz, '%.8f') '];']);
disp(' ');


% perform calibration
if (true)
    %temperature calibration (using filtered temperature to minimize noise)
    T_padded = [data.SparkFun6DOF.temperature(1)*ones(250,1); data.SparkFun6DOF.temperature; ...
        data.SparkFun6DOF.temperature(end)*ones(250,1)];
    T = filter2(ones(501,1)/501, T_padded, 'valid');
    acc1_x = (data.SparkFun6DOF.acc_x - (a1_b_Ax*T + a1_b_Bx)) ./ (a1_g_Ax*T + a1_g_Bx);
    acc1_y = (data.SparkFun6DOF.acc_y - (a1_b_Ay*T + a1_b_By)) ./ (a1_g_Ay*T + a1_g_By);
    acc1_z = (data.SparkFun6DOF.acc_z - (a1_b_Az*T + a1_b_Bz)) ./ (a1_g_Az*T + a1_g_Bz);
    
    w1_x = (data.SparkFun6DOF.w_x - (w1_b_Ax*T + w1_b_Bx)) ./ (w1_g_Ax*T + w1_g_Bx);
    w1_y = (data.SparkFun6DOF.w_y - (w1_b_Ay*T + w1_b_By)) ./ (w1_g_Ay*T + w1_g_By);
    w1_z = (data.SparkFun6DOF.w_z - (w1_b_Az*T + w1_b_Bz)) ./ (w1_g_Az*T + w1_g_Bz);
    
    if (false)
        % Unit test
        acc1_x_ref = (data.SparkFun6DOF.acc_x - B1(1)) / G1(1);
        acc1_y_ref = (data.SparkFun6DOF.acc_y - B1(2)) / G1(2);
        acc1_z_ref = (data.SparkFun6DOF.acc_z - B1(3)) / G1(3);

        w1_x_ref = (data.SparkFun6DOF.w_x - w1_Bx) / w1_Gx;
        w1_y_ref = (data.SparkFun6DOF.w_y - w1_By) / w1_Gy;
        w1_z_ref = (data.SparkFun6DOF.w_z - w1_Bz) / w1_Gz;    

        notNan = ~isnan(acc1_x);
        ax_error = mean(acc1_x(notNan) - acc1_x_ref(notNan))
        ay_error = mean(acc1_y(notNan) - acc1_y_ref(notNan))
        az_error = mean(acc1_z(notNan) - acc1_z_ref(notNan))

        wx_error = mean(w1_x(notNan) - w1_x_ref(notNan))
        wy_error = mean(w1_y(notNan) - w1_y_ref(notNan))
        wz_error = mean(w1_z(notNan) - w1_z_ref(notNan))
    end
else
    %simple calibration
    acc1_x = (data.SparkFun6DOF.acc_x - B1(1)) / G1(1);
    acc1_y = (data.SparkFun6DOF.acc_y - B1(2)) / G1(2);
    acc1_z = (data.SparkFun6DOF.acc_z - B1(3)) / G1(3);
    
    w1_x = (data.SparkFun6DOF.w_x - w1_Bx) / w1_Gx;
    w1_y = (data.SparkFun6DOF.w_y - w1_By) / w1_Gy;
    w1_z = (data.SparkFun6DOF.w_z - w1_Bz) / w1_Gz;
end

acc2_x = (data.InertiaLink.acc_x - B2(1)) / G2(1);
acc2_y = (data.InertiaLink.acc_y - B2(2)) / G2(2);
acc2_z = (data.InertiaLink.acc_z - B2(3)) / G2(3);

w2_x = (data.InertiaLink.w_x - w2_Bx) / w2_Gx;
w2_y = (data.InertiaLink.w_y - w2_By) / w2_Gy;
w2_z = (data.InertiaLink.w_z - w2_Bz) / w2_Gz;


if (false) %debugplotting for calibration

    figure(1);
    clf;
    subplot(3,1,1);
    plot(SFtime, acc1_x, 'c', ILtime, acc2_x, 'b--', 'linewidth', 1); 
    hold on;
    minY = min([min(acc1_x) min(acc2_x)]);
    maxY = max([max(acc1_x) max(acc2_x)]);
    for i = 1:size(calibTimes,2)
        plot([calibTimes(1,i) calibTimes(1,i)], [minY maxY], 'k:');
        plot([calibTimes(2,i) calibTimes(2,i)], [minY maxY], 'k:');
        
        plot([calibTimes2(1,i) calibTimes2(1,i)], [minY maxY], 'k:');
        plot([calibTimes2(2,i) calibTimes2(2,i)], [minY maxY], 'k:');        
    end
    for i = 1:size(angleCalibTimes,2)
        plot([angleCalibTimes(1,i) angleCalibTimes(1,i)], [minY maxY], 'r:');
        plot([angleCalibTimes(end,i) angleCalibTimes(end,i)], [minY maxY], 'r:');
        
        plot([angleCalibTimes2(1,i) angleCalibTimes2(1,i)], [minY maxY], 'r:');
        plot([angleCalibTimes2(end,i) angleCalibTimes2(end,i)], [minY maxY], 'r:');        
    end    
    axis([minCalibTime, maxCalibTime, minY, maxY])
    ylabel('x');

    subplot(3,1,2);
    plot(SFtime, acc1_y, 'c', ILtime, acc2_y, 'b--', 'linewidth', 1); 
    hold on;
    minY = min([min(acc1_y) min(acc2_y)]);
    maxY = max([max(acc1_y) max(acc2_y)]);
    for i = 1:size(calibTimes,2)
        plot([calibTimes(1,i) calibTimes(1,i)], [minY maxY], 'k:');
        plot([calibTimes(2,i) calibTimes(2,i)], [minY maxY], 'k:');
        
        plot([calibTimes2(1,i) calibTimes2(1,i)], [minY maxY], 'k:');
        plot([calibTimes2(2,i) calibTimes2(2,i)], [minY maxY], 'k:');          
    end
    for i = 1:size(angleCalibTimes,2)
        plot([angleCalibTimes(1,i) angleCalibTimes(1,i)], [minY maxY], 'r:');
        plot([angleCalibTimes(end,i) angleCalibTimes(end,i)], [minY maxY], 'r:');
                
        plot([angleCalibTimes2(1,i) angleCalibTimes2(1,i)], [minY maxY], 'r:');
        plot([angleCalibTimes2(end,i) angleCalibTimes2(end,i)], [minY maxY], 'r:');   
    end    
    axis([minCalibTime, maxCalibTime, minY, maxY])
    ylabel('y');

    subplot(3,1,3);
    plot(SFtime, acc1_z, 'c', ILtime, acc2_z, 'b--', 'linewidth', 1); 
    hold on;
    minY = min([min(acc1_z) min(acc2_z)]);
    maxY = max([max(acc1_z) max(acc2_z)]);
    for i = 1:size(calibTimes,2)
        plot([calibTimes(1,i) calibTimes(1,i)], [minY maxY], 'k:');
        plot([calibTimes(2,i) calibTimes(2,i)], [minY maxY], 'k:');
        
        plot([calibTimes2(1,i) calibTimes2(1,i)], [minY maxY], 'k:');
        plot([calibTimes2(2,i) calibTimes2(2,i)], [minY maxY], 'k:');            
    end
    for i = 1:size(angleCalibTimes,2)
        plot([angleCalibTimes(1,i) angleCalibTimes(1,i)], [minY maxY], 'r:');
        plot([angleCalibTimes(end,i) angleCalibTimes(end,i)], [minY maxY], 'r:');
        
        plot([angleCalibTimes2(1,i) angleCalibTimes2(1,i)], [minY maxY], 'r:');
        plot([angleCalibTimes2(end,i) angleCalibTimes2(end,i)], [minY maxY], 'r:');   
    end    
    ylabel('z');
    axis([minCalibTime, maxCalibTime, minY, maxY])
    pause(0.2);
    
    
    figure(2);
    clf;
    subplot(3,1,1);
    plot(Ktime, gyro(:,1), 'k', SFtime-SFdelay, w1_x, 'c', ILtime-ILdelay, w2_x, 'b--', 'linewidth', 1); 
    hold on;
    minY = min([min(w1_x) min(w2_x)]);
    maxY = max([max(w1_x) max(w2_x)]);
    for i = 1:size(calibTimes,2)
        plot([calibTimes(1,i) calibTimes(1,i)], [minY maxY], 'k:');
        plot([calibTimes(2,i) calibTimes(2,i)], [minY maxY], 'k:');
    end
    for i = 1:size(angleCalibTimes,2)
        plot([angleCalibTimes(1,i) angleCalibTimes(1,i)], [minY maxY], 'r:');
        plot([angleCalibTimes(end,i) angleCalibTimes(end,i)], [minY maxY], 'r:');
    end    
    axis([minTime, maxTime, minY, maxY])
    ylabel('x');

    subplot(3,1,2);
    plot(Ktime, gyro(:,2), 'k', SFtime-SFdelay, w1_y, 'c', ILtime-ILdelay, w2_y, 'b--', 'linewidth', 1); 
    hold on;
    minY = min([min(w1_y) min(w2_y)]);
    maxY = max([max(w1_y) max(w2_y)]);
    for i = 1:size(calibTimes,2)
        plot([calibTimes(1,i) calibTimes(1,i)], [minY maxY], 'k:');
        plot([calibTimes(2,i) calibTimes(2,i)], [minY maxY], 'k:');
    end
    for i = 1:size(angleCalibTimes,2)
        plot([angleCalibTimes(1,i) angleCalibTimes(1,i)], [minY maxY], 'r:');
        plot([angleCalibTimes(end,i) angleCalibTimes(end,i)], [minY maxY], 'r:');
    end    
    axis([minTime, maxTime, minY, maxY])
    ylabel('y');

    subplot(3,1,3);
    plot(Ktime, gyro(:,3), 'k', SFtime-SFdelay, w1_z, 'c', ILtime-ILdelay, w2_z, 'b--', 'linewidth', 1); 
    hold on;
    minY = min([min(w1_z) min(w2_z)]);
    maxY = max([max(w1_z) max(w2_z)]);
    for i = 1:size(calibTimes,2)
        plot([calibTimes(1,i) calibTimes(1,i)], [minY maxY], 'k:');
        plot([calibTimes(2,i) calibTimes(2,i)], [minY maxY], 'k:');
    end
    for i = 1:size(angleCalibTimes,2)
        plot([angleCalibTimes(1,i) angleCalibTimes(1,i)], [minY maxY], 'r:');
        plot([angleCalibTimes(end,i) angleCalibTimes(end,i)], [minY maxY], 'r:');
    end    
    ylabel('z');
    axis([minTime, maxTime, minY, maxY])
    pause(0.2);

    if (false) %test calibration parameters (run calibration with corrected data)
        w_x_error_SF = w1_x(SFrange) - w_x_in_SF;
        w_y_error_SF = w1_y(SFrange) - w_y_in_SF;
        w_z_error_SF = w1_z(SFrange) - w_z_in_SF;

        w_x_error_IL = w2_x(ILrange) - w_x_in_IL;
        w_y_error_IL = w2_y(ILrange) - w_y_in_IL;
        w_z_error_IL = w2_z(ILrange) - w_z_in_IL;

        [w1_Gx, w1_Bx] = fitRegressionLine(w_x_in_SF(x_filt_SF), w_x_in_SF(x_filt_SF)+w_x_error_SF(x_filt_SF));
        [w2_Gx, w2_Bx] = fitRegressionLine(w_x_in_IL(x_filt_IL), w_x_in_IL(x_filt_IL)+w_x_error_IL(x_filt_IL));
        [w1_Gy, w1_By] = fitRegressionLine(w_y_in_SF(y_filt_SF), w_y_in_SF(y_filt_SF)+w_y_error_SF(y_filt_SF));
        [w2_Gy, w2_By] = fitRegressionLine(w_y_in_IL(y_filt_IL), w_y_in_IL(y_filt_IL)+w_y_error_IL(y_filt_IL));
        [w1_Gz, w1_Bz] = fitRegressionLine(w_z_in_SF(z_filt_SF), w_z_in_SF(z_filt_SF)+w_z_error_SF(z_filt_SF));
        [w2_Gz, w2_Bz] = fitRegressionLine(w_z_in_IL(z_filt_IL), w_z_in_IL(z_filt_IL)+w_z_error_IL(z_filt_IL));    
    end
    
    figure(3);
    clf;
    subplot(3,1,1);
    minY = -pi;
    maxY = pi;
    plot(w_x_in_SF(x_filt_SF), w_x_error_SF(x_filt_SF), 'b.', w_x_in_IL(x_filt_IL), w_x_error_IL(x_filt_IL), 'm.', 'linewidth', 1, 'markersize', 3); 
    hold on;
    x_vals = [min(w_x_in_SF(x_filt_SF)), max(w_x_in_SF(x_filt_SF))];
    plot(x_vals, (w1_Gx-1)*x_vals + w1_Bx, 'b', x_vals, (w2_Gx-1)*x_vals + w2_Bx, 'm', 'linewidth', 1);
    ylabel('w_x error (rad/s)');
    xlabel('reference w_x (rad/s)');
    
    
    subplot(3,1,2);
    plot(w_y_in_SF(y_filt_SF), w_y_error_SF(y_filt_SF), 'b.', w_y_in_IL(y_filt_IL), w_y_error_IL(y_filt_IL), 'm.', 'linewidth', 1, 'markersize', 3); 
    hold on;
    y_vals = [min(w_y_in_SF(y_filt_SF)), max(w_y_in_SF(y_filt_SF))];
    plot(y_vals, (w1_Gy-1)*y_vals + w1_By, 'b', y_vals, (w2_Gy-1)*y_vals + w2_By, 'm', 'linewidth', 1);
    ylabel('w_y error (rad/s)');
    xlabel('reference w_y (rad/s)');
    %ylabel('w_y (rad/s)');

    subplot(3,1,3);
    plot(w_z_in_SF(z_filt_SF), w_z_error_SF(z_filt_SF), 'b.', w_z_in_IL(z_filt_IL), w_z_error_IL(z_filt_IL), 'm.', 'linewidth', 1, 'markersize', 3); 
    hold on;
    z_vals = [min(w_z_in_SF(z_filt_SF)), max(w_z_in_SF(z_filt_SF))];
    plot(z_vals, (w1_Gz-1)*z_vals + w1_Bz, 'b', z_vals, (w2_Gz-1)*z_vals + w2_Bz, 'm', 'linewidth', 1);
    ylabel('w_z error (rad/s)');
    xlabel('reference w_z (rad/s)');
    %ylabel('w_z (rad/s)');

    pause(0.2);
end
