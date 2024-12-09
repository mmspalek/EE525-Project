% Nolan Egging
% EE525 - Draft 3
% Data Analysis
% Due Nov 22

close all;
clear;
clc;

% system variables
T = 0.02;
horizon = 38; 
t = 0:T:horizon;

% corrective factor
cf = 1.5;

% state space model
A_ct = [0 1 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 1 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 1;
        0 0 0 0 0 0];
B_ct = [0 0 0;
        1 0 0;
        0 0 0;
        0 1 0;
        0 0 0;
        0 0 1];
C = [1 0 1 0 0 0;
     1 0 0 0 1 0];
D = [0 0 0;
     0 0 0];

% CT -> DT
trick1 = expm([A_ct, B_ct; zeros(3,9)]*T); % found on wikipedia
A = trick1(1:6, 1:6);
B = trick1(1:6, 7:9);

% raw data
mpu_table = readtable('MPU6050_upDn.xlsx');
mtb_table = readtable('MATLAB_upDn.xlsx');
mpu_raw = rmmissing(mpu_table{:, 1:15}) * (9.81 / 2^14);
mtb_raw = rmmissing(mtb_table{:, 1:15});

acc1 = mpu_raw(44:end, 4)';
acc2 = mtb_raw(30:end, 4)';

acc1(1:1000) = acc1(1:1000) - mean(acc1(1:2/T));
acc1(1001:end) = acc1(1001:end) - mean(acc1((33/T):(35/T)));
acc1 = acc1(1:length(t));
acc2 = acc2 - mean(acc2(1:2/T));
acc2(1:1000) = acc2(1:1000) - mean(acc2(1:2/T));
acc2(1001:end) = acc2(1001:end) - mean(acc2((33/T):(35/T)));
acc2 = acc2(1:length(t));

% plot raw data
figure(1);
plot(t, acc1);
hold on;
plot(t, acc2);
hold off;
legend("MPU Acc.", "MATLAB Acc.");

% integrate data
single_int = zeros(2, length(t));
y = zeros(2, length(t));
for i = 2:1:length(t)
    single_int(1, i) = single_int(1, i-1) + acc1(i-1)*T;
    single_int(2, i) = single_int(2, i-1) + acc2(i-1)*T;
    y(1, i) = y(1, i-1) + single_int(1, i-1)*T;
    y(2, i) = y(2, i-1) + single_int(2, i-1)*T;
end

% constuct ideal input model
[a, p] = getModel(T, horizon);
u(1, :) = a;

% apply penrose pseudo inverse
xls = zeros(6, length(t));
k = 3;
Ok = [C; C*A; C*A*A];
Tk = [D, zeros(size(D)), zeros(size(D));
      C*B, D zeros(size(D));
      C*A*B, C*B, D];
for j = 1:1:length(t)
    if j+k-1 == length(t)
        break;
    end
    uk = [u(1, j:j+k-1); zeros(1, k); zeros(1, k)];
    uk = uk(:);
    yk = y(:, j:j+k-1);
    yk = yk(:);
    xls(:, j) = pinv(Ok) * (yk - (Tk*uk));
end

% plot height measurement data
figure(2);
plot(t, y(1, :));
hold on;
plot(t, y(2, :));
plot(t, xls(1, :));
hold off;
legend("y1", "y2", "xls");

% put entries into a compressed data format
rd1 = zeros(15, length(t));
rd2 = zeros(15, length(t));
mpu_shifts = [22, 81, 62, 44, 56, 82, 65, 49, 39, 53, 37, 77, 69, 72, 43];
mtb_shifts = [1, 53, 44, 30, 59, 41, 48, 70, 54, 42, 48, 29, 61, 60, 58];
for i = 1:1:15
    rd1(i, :) = mpu_raw(mpu_shifts(i):mpu_shifts(i)+length(t)-1, i);
    rd2(i, :) = mtb_raw(mtb_shifts(i):mtb_shifts(i)+length(t)-1, i);
end

% processes each run to find average
sims = 15;
sum_y1 = zeros(1, length(t));
sum_y2 = zeros(1, length(t));
sum_x1_ls = zeros(1, length(t));
for i = 1:1:sims
    acc1 = rd1(i, :);
    acc2 = rd2(i, :);

    % remove mean
    acc1(1:1000) = acc1(1:1000) - mean(acc1(1:2/T));
    acc1(1001:end) = acc1(1001:end) - mean(acc1((33/T):(35/T)));
    acc1 = acc1(1:length(t));
    acc2 = acc2 - mean(acc2(1:2/T));
    acc2(1:1000) = acc2(1:1000) - mean(acc2(1:2/T));
    acc2(1001:end) = acc2(1001:end) - mean(acc2((33/T):(35/T)));
    acc2 = acc2(1:length(t));

    % integrate data
    single_int = zeros(2, length(t));
    y = zeros(2, length(t));
    for j = 2:1:length(t)
        single_int(1, j) = single_int(1, j-1) + acc1(j-1)*T;
        single_int(2, j) = single_int(2, j-1) + acc2(j-1)*T;
        y(1, j) = y(1, j-1) + single_int(1, j-1)*T;
        y(2, j) = y(2, j-1) + single_int(2, j-1)*T;
    end

    % apply LS estimate
    xls = zeros(6, length(t));
    k = 3;
    for j = 1:1:length(t)
        if j+k-1 == length(t)
            break;
        end
        uk = [u(1, j:j+k-1); zeros(1, k); zeros(1, k)];
        uk = uk(:);
        yk = y(:, j:j+k-1);
        yk = yk(:);
        xls(:, j) = pinv(Ok) * (yk - (Tk*uk));
    end

    sum_y1 = sum_y1 + y(1, :);
    sum_y2 = sum_y2 + y(2, :);
    sum_x1_ls = sum_x1_ls + xls(1, :);

end

% calculate and plot averages
avg_y1 = sum_y1 / sims;
avg_y2 = sum_y2 / sims;
avg_x1_ls = sum_x1_ls / sims;
figure(3);
plot(t, avg_y1, t, avg_y2, t, 1.5*avg_x1_ls);
legend("Average y1", "Average y2", "Average LS x1");