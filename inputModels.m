% Nolan Egging
% EE525 - Draft 3
% Input Models
% Due Nov 22

close all;
clear;
clc;

% system variables
T = 0.02;
t = (0:T:40);

% read-in data and get average
mpu_table = readtable('MPU6050_upDn.xlsx');
mtb_table = readtable('MATLAB_upDn.xlsx');
mpu_raw = rmmissing(mpu_table{:, 1:15}) * (9.81 / 2^14);
mtb_raw = rmmissing(mtb_table{:, 1:15});
avg = (mean(mpu_raw, 2) + mean(mtb_raw, 2)) * 0.5;
%figure(1);
%hold on;
%plot(mpu_raw(22:end, 1));
%plot(mpu_raw(81:end, 2));
%plot(mpu_raw(62:end, 3));
%plot(mpu_raw(44:end, 4));
%plot(mpu_raw(56:end, 5));
%plot(mpu_raw(82:end, 6));
%plot(mpu_raw(65:end, 7));
%plot(mpu_raw(49:end, 8));
%plot(mpu_raw(39:end, 9));
%plot(mpu_raw(53:end, 10));
%plot(mpu_raw(37:end, 11));
%plot(mpu_raw(77:end, 12));
%plot(mpu_raw(69:end, 13));
%plot(mpu_raw(72:end, 14));
%plot(mpu_raw(43:end, 15));
%plot(mtb_raw(1:end, 1));
%plot(mtb_raw(53:end, 2));
%plot(mtb_raw(44:end, 3));
%plot(mtb_raw(30:end, 4));
%plot(mtb_raw(59:end, 5));
%plot(mtb_raw(41:end, 6));
%plot(mtb_raw(48:end, 7));
%plot(mtb_raw(70:end, 8));
%plot(mtb_raw(54:end, 9));
%plot(mtb_raw(42:end, 10));
%plot(mtb_raw(48:end, 11));
%plot(mtb_raw(29:end, 12));
%plot(mtb_raw(61:end, 13));
%plot(mtb_raw(60:end, 14));
%plot(mtb_raw(58:end, 15));
%hold off;

% clean up data
rmavg_data = zeros(1, length(avg));
rmavg_data(1:1000) = avg(1:1000) - mean(avg(1:2/T)); % remove mean of first 2 seconds
rmavg_data(1001:end) = avg(1001:end) - mean(avg((34.5/T):(36.5/T))); % remove mean of last 2 seconds
clean_data = rmavg_data(1:1:length(t));

% create model
[a_model, p_model] = getModel(T, (length(clean_data)-1)*T);

% plot model against mpu
figure(2);
plot(t, clean_data);
hold on;
plot(t, a_model(1:length(t)));
hold off;
lgd1 = legend("30 Run Average", "Model");
lgd1.Location = "southeast";
title("Sensor Acceleration versus Model");
xlabel("Time (s)");
ylabel("Input Acceleration (m/s/s)");
xlim([0 40]);

% integrate model
model_int = (0:T:(length(clean_data)-1)*T);
model_double = (0:T:(length(clean_data)-1)*T);
for i = 2:1:length(model_int)
    model_int(i) = model_int(i-1) + a_model(i-1)*T;
    model_double(i) = model_double(i-1) + model_int(i-1)*T;
end
% plot model integral and double integral
figure(3);
subplot(211);
plot(t, model_int);
title("Integral of Model: Model Velocity");
xlabel("Time (s)");
ylabel("Velocity (m/s)");
xlim([0 39]);
ylim([-0.6 0.6]);
subplot(212);
plot(t, model_double);
title("Double Integral of Model: Model Position");
xlabel("Time (s)");
ylabel("Position (m)");
xlim([0 39]);
ylim([-0.5 4.0]);

% plot p_model and double integral
figure(4);
plot(t, model_double, t, p_model);
title("Double Integral of Acceleration Model vs Position Model");
xlabel("Time (s)");
ylabel("Position (m)");
lgd2 = legend("Double Integral of Acc. Model", "Position Model");
lgd2.Location = "east";