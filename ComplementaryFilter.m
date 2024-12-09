% Nolan Egging
% EE525 - Final
% Complementary
% Due Dec 12th

close all; clear; clc;

% system variables
dt = 0.02;
horizon = 38; 
t = 0:dt:horizon;
steps = length(t);
epsilon = 1e-9;

% read-in data and get average
mpu_table = readtable('MPU6050_upDn.xlsx');
mtb_table = readtable('MATLAB_upDn.xlsx');
mpu_raw = rmmissing(mpu_table{:, 1:15}) * (9.81 / 2^14);
mtb_raw = rmmissing(mtb_table{:, 1:15});

% clean up data
%acc1 = mpu_raw(22:end, 1)';
%acc2 = mtb_raw(1:end, 1)';
%acc1 = mpu_raw(81:end, 2)';
%acc2 = mtb_raw(53:end, 2)';
acc1 = mpu_raw(62:end, 3)';
acc2 = mtb_raw(44:end, 3)';
%acc1 = mpu_raw(44:end, 4)';
%acc2 = mtb_raw(30:end, 4)';
%acc1 = mpu_raw(56:end, 5)';
%acc2 = mtb_raw(59:end, 5)';
%acc1 = mpu_raw(82:end, 6)';
%acc2 = mtb_raw(41:end, 6)';
acc1(1:1000) = acc1(1:1000) - mean(acc1(1:2/dt));
acc1(1001:end) = acc1(1001:end) - mean(acc1((33/dt):(35/dt)));
acc1 = acc1(1:length(t));
acc2 = acc2 - mean(acc2(1:2/dt));
acc2(1:1000) = acc2(1:1000) - mean(acc2(1:2/dt));
acc2(1001:end) = acc2(1001:end) - mean(acc2((33/dt):(35/dt)));
acc2 = acc2(1:length(t));
z =[acc2; acc1];

% characterized noises
mpu_var = 4.4e-4;
mtb_wn_var = 2.69e-5;
mtb_wn_beta = 100;
mtb_gm_var = 5.5e-4;
mtb_gm_beta = 0.3;
unity_wn = 1;

% loop through process variables
proc_vars = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
proc_beta = 1/(5*dt);
error = zeros(6, steps);
for trial = 1:1:6
    
    % set proccess variation
    proc_var = proc_vars(trial);

    % system variables
    n = 3;
    A = eye(n) .* [-mtb_wn_beta, -mtb_gm_beta, -proc_beta];
    B = eye(n) .* [sqrt(2*mtb_wn_var), sqrt(2*mtb_gm_var), sqrt(2*proc_var)];
    H = [1 1 1;
         0 0 1];
    W = eye(n) .* unity_wn;
    
    % van loan discretization
    L = [-A, B*W*B'; zeros(n), A'] * dt;
    E = expm(L);
    Phi = E(n+1:end, n+1:end)';
    Q = Phi * E(1:n, n+1:end);
    
    % sensor covariance
    R = eye(2) .* [epsilon, mpu_var];
    
    % initial guesses
    xhatminus = zeros(n, 1);
    Pminus = eye(3) .* [mtb_wn_var, mtb_gm_var, proc_var];
    
    % estimate and error placeholders
    xhat = zeros(n, steps);
    P = zeros(n, n, steps);
    
    % Kalman filter loop
    for i = 1:1:steps
    
        % step 1: compute Kalman gain
        K = Pminus * H' / (H*Pminus*H'+ R);
    
        % step 2: get measurement and update estimate
        xhat(:, i) = xhatminus + K*(z(:, i) - H*xhatminus);
    
        % step 3: update error covariance
        P(:, :, i) = (eye(n) - K*H) * Pminus;
    
        % step 4: project ahead
        xhatminus = (Phi * xhat(:, i));
        Pminus = Phi*P(:, :, i)*Phi' + Q;
    
    end
    est_acc = xhat(3, :);

    error(trial, :) = P(3, 3, :);

end

% plot acceleration data
figure(1);
plot(t, acc1, t, acc2, t, est_acc);
legend("MPU Acc.", "MATLAB Acc.", "Comp. Est");
xlabel("Time (s)");
ylabel("Acceleration (m/s^2) - Gravity");
xlim([0 38]);

% integrate data
single_int = zeros(3, length(t));
pos = zeros(3, length(t));
for i = 2:1:steps
    single_int(1, i) = single_int(1, i-1) + acc1(i-1)*dt;
    single_int(2, i) = single_int(2, i-1) + acc2(i-1)*dt;
    single_int(3, i) = single_int(3, i-1) + est_acc(i-1)*dt;
    pos(1, i) = pos(1, i-1) + single_int(1, i-1)*dt;
    pos(2, i) = pos(2, i-1) + single_int(2, i-1)*dt;
    pos(3, i) = pos(3, i-1) + single_int(3, i-1)*dt;
end

% plot double integral data
figure(3);
plot(t, pos(1, :), t, pos(2, :), t, pos(3, :));
legend("MPU", "MATLAB", "Estimate");
xlabel("Time (s)");
ylabel("Position (m)");
xlim([0 38]);

% plots convergence of error
figure(4);
plot(t, error(1, :), t, error(2, :), t, error(3, :), t, error(4, :), t, error(5, :), t, error(6, :));
legend("\sigma^2 = 10^{-5}", "\sigma^2 = 10^{-4}", "\sigma^2 = 10^{-3}", "\sigma^2 = 10^{-2}", "\sigma^2 = 10^{-1}", "\sigma^2 = 1");
ylabel('Acceleration MSE ((m/s^2)^2)'); xlabel("Time (s)");
xlim([0 35]);


% integrate error covariance
single_int_var = zeros(1, length(t));
pos_var = zeros(1, length(t));
for i = 2:1:steps
    single_int_var(1, i) = single_int_var(1, i-1) + error(6, i-1)*dt;
    pos_var(1, i) = pos_var(1, i-1) + single_int_var(1, i-1)*dt;
end

% plots filter error covariance of position
figure(5);
plot(t, pos_var);
ylabel('Position MSE (m^2)'); xlabel("Time (s)");

% display heights and drifts
t1 = 11.1;
fprintf("%.2f \\\\ \\hline \n", pos(3, t1/dt));
t1 = 32.4;
fprintf("%.2f \\\\ \\hline \n", pos(3, t1/dt));