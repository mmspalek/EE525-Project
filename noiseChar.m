% Nolan Egging
% EE525 - Draft 2
% Noise Characterization
% Due Nov 1
%
% Referenced Prof. McKell's periodogram_ext.m

close all;
clear;
clc;

% 1 -> MPU6050
% 2 -> MATLAB Sensor

% get data
table1 = readtable('mpu6050_data.xlsx');
table2 = readtable('matlab_data3.xlsx');
raw_data1 = rmmissing(table1{:, 2:11});
raw_data2 = rmmissing(table2{:, 2:11});
X1 = raw_data1 * (9.81/(2^14)); % convert to acceleration
X2 = raw_data2; % matlab outputs in m/s^2
dt = 0.02; % fixed sample rate for both MPU6050 and MATLAB
[samples1, runs1] = size(X1);
[samples2, runs2] = size(X2);
t1 = dt:dt:samples1*dt;
t2 = dt:dt:samples2*dt;

% plot raw data and average raw data
figure(1);
subplot(211);
plot(t1', X1(:, 10));
hold on;
plot(t2', X2(:, 10));
hold off;
title("Raw Signals for Last Run");
xlabel("Time (s)");
xlim([0 60]);
ylim([9.7 10.0]);
ylabel("Acceleration + Gravity (m/s^2)");
legend("MPU6050", "MATLAB");
subplot(212);
plot(t1', mean(X1, 2));
hold on;
plot(t2', mean(X2, 2));
hold off;
title("Mean Raw Signal");
xlabel("Time (s)");
xlim([0 60]);
ylim([9.7 10.0]);
ylabel("Acceleration + Gravity (m/s^2)");
X1_mean = mean(X1);
X2_mean = mean(X2);
X1 = X1 - X1_mean;
X2 = X2 - X2_mean;
disp("Mean Gravity (MPU6050) = " + mean(X1_mean));
disp("Mean Gravity (MATLAB) = " + mean(X2_mean));

% autocorrelation variables
T1 = samples1*dt;
T2 = samples2*dt;
max_tau = 40;
tau = (0:dt:max_tau);

% calculate autocorrelation
Vx1 = zeros(runs1, length(tau));
for run = 1:1:runs1
    for i = 1:1:length(tau)
        Vx1(run, i) = (1/(T1 - tau(i))) * (X1(1:samples1-i+1, run)' * X1(i:samples1, run) * dt);
    end
end
Vx2 = zeros(runs2, length(tau));
for run = 1:1:runs2
    for i = 1:1:length(tau)
        Vx2(run, i) = (1/(T2 - tau(i))) * (X2(1:samples2-i+1, run)' * X2(i:samples2, run) * dt);
    end
end

% plot autocorrelation
figure(2);
subplot(211);
plot(tau, Vx1(10, :));
title("Autocorrelation of Last Run (MPU6050)");
xlabel("\tau (s)");
ylabel("Acceleration^2 ((m/s^2)^2)");
xlim([-1 max_tau+1]);
ylim([-0.1*Vx1(10, 1) 1.1*Vx1(10, 1)]);
text(2, Vx1(10, 1)*0.8, "Sample Variance = " + Vx1(10, 1));
subplot(212);
plot(tau, Vx2(10, :), "Color", "#D95319");
title("Autocorrelation of Last Run (MATLAB)");
xlabel("\tau (s)");
ylabel("Acceleration^2 ((m/s^2)^2)");
xlim([-1 max_tau+1]);
ylim([-0.1*Vx2(10, 1) 1.1*Vx2(10, 1)]);
text(2, Vx2(10, 1)*0.8, "Sample Variance = " + Vx2(10, 1));

% plot mean autocorrelation
figure(3);
subplot(211);
plot(tau, mean(Vx1));
title("Mean Autocorrelation (MPU6050)");
xlabel("\tau (s)");
ylabel("Acceleration^2 ((m/s^2)^2)");
xlim([-1 max_tau+1]);
mean_var1 = mean(Vx1(:, 1));
text(2, mean_var1*0.8, "Mean sample Variance = " + mean_var1);
ylim([-0.1*mean_var1 1.1*mean_var1]);
disp("Sample Variance from Mean Autocorrelation = " + mean_var1);
subplot(212);
plot(tau, mean(Vx2), "Color", "#D95319");
title("Mean Autocorrelation (MATLAB)");
xlabel("\tau (s)");
ylabel("Acceleration^2 ((m/s^2)^2)");
xlim([-1 max_tau+1]);
mean_var2 = mean(Vx2(:, 1));
text(2, mean_var2*0.8, "Mean sample Variance = " + mean_var2);
ylim([-0.1*mean_var2 1.1*mean_var2]);
disp("Sample Variance from Mean Autocorrelation = " + mean_var2);

% create a periodogram
figure(4);
[psd1, f1] = periodogram(X1(:, 10), hamming(samples1), samples1, 1/dt);
[psd2, f2] = periodogram(X2(:, 10), hamming(samples2), samples2, 1/dt);
% Units of periodogram are (sampled units)^2/Hz
% so (m/s/s)^2/Hz
subplot(211);
plot(f1, db(psd1,'power'));
hold on;
plot(f2, db(psd2,'power'));
hold off;
ylim([-100 0]);
title("Power Spectral Density for Last Run");
xlabel('Analog frequency [Hz]');
ylabel('PSD estimate [dB_{m/s^2}/Hz]');
legend("MPU6050", "MATLAB");

% create average periodogram
psd1_mean = zeros(length(f1), 1);
psd2_mean = zeros(length(f2), 1);
for i = 1:1:runs1
    [psd1_inc,] = periodogram(X1(:, i), hamming(samples1), samples1, 1/dt);
    [psd2_inc,] = periodogram(X2(:, i), hamming(samples2), samples2, 1/dt);
    psd1_mean = psd1_mean + psd1_inc;
    psd2_mean = psd2_mean + psd2_inc;
end
psd1_mean = psd1_mean / runs1;
psd2_mean = psd2_mean / runs2;
subplot(212);
plot(f1, db(psd1_mean,'power'), "Color", "#0072BD");
hold on;
plot(f2, db(psd2_mean,'power'), "Color", "#D95319");
hold off;
ylim([-100 0]);
title("Mean Power Spectral Density");
xlabel('Analog frequency [Hz]');
ylabel('PSD estimate [dB_{m/s^2}/Hz]');
disp("Mean PSD (MPU6050) = " + mean(db(psd1_mean,'power')));
disp("Mean PSD (MATLAB0) = " + mean(db(psd2_mean,'power')));

% integrates noise over time for process
figure(5);
X1_cum = zeros(samples2, 10);
X2_cum = zeros(samples2, 10);
for i = 1:1:runs1
    X1_cum(:, i) = cumtrapz(cumtrapz(X1(:, i))) * dt * dt;
end
for i = 1:1:runs2
    X2_cum(:, i) = cumtrapz(cumtrapz(X2(:, i))) * dt * dt;
end
subplot(211); 
plot(t1, X1_cum(:, 10));
hold on;
plot(t2, X2_cum(:, 10));
hold off;
title("Double Integral for Last Run");
xlabel("Time (s)");
xlim([0 60]);
ylabel("Position (m)");
legend("MPU6050", "MATLAB");

% plots mean integrated noise
X1_mean_cum = mean(X1_cum, 2);
X2_mean_cum = mean(X2_cum, 2);
subplot(212); 
plot(t1, X1_mean_cum);
hold on;
plot(t2, X2_mean_cum);
hold off;
title("Mean of Double Integral");
xlabel("Time (s)");
xlim([0 60]);
ylabel("Position (m)");

% find variance over time
figure(6);
X1_var_cum = mean(X1_cum.^2, 2);
X2_var_cum = mean(X2_cum.^2, 2);
plot(t1, X1_var_cum);
hold on;
plot(t2, X2_var_cum);
hold off;
title("Variance of Double Integral");
xlabel("Time (s)");
xlim([0 60]);
ylabel("Position^2 (m^2)");
legend("MPU6050", "MATLAB");

% plot fit autocorrelation for the MATLAB sensor
beta_wn = 25; % 100
beta_slow = 0.3; % 0.3
var_wn = 2.69e-5;
var_slow = 0.55e-5; % 0.55e-5
Rx = (var_wn*exp(-beta_wn*tau)) + (var_slow*exp(-beta_slow*tau));
figure(7);
plot(tau, mean(Vx2), "Color", "#D95319");
hold on;
plot(tau, Rx, 'k');
hold off;
title("Mean Autocorrelation versus Best Fit (MATLAB)");
xlabel("\tau (s)");
ylabel("Acceleration^2 ((m/s^2)^2)");
legend("Experimental", "Fit");
text(2, mean_var2*0.8, "Mean sample Variance = " + mean_var2);
xlim([-1 41]);
ylim([-0.1*mean_var2 1.1*mean_var2]);

% plots moving mean
figure(8)
plot(t1', movmean(mean(X1, 2), 300));
hold on;
plot(t2', movmean(mean(X2, 2), 300));
hold off;
title("Moving Average of 300 Samples of Raw Data");
xlabel("Time (s)");
xlim([0 60]);
ylabel("Acceration (m/s^2)");
legend("MPU6050", "MATLAB");
