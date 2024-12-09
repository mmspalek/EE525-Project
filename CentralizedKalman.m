% Nolan Egging
% EE525 - Final
% Centralized Kalman
% Due Dec 12th

close all; clear; clc;

% system variables
dt = 0.02;
horizon = 38; 
t = 0:dt:horizon;
steps = length(t);

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
%acc1 = mpu_raw(62:end, 3)';
%acc2 = mtb_raw(44:end, 3)';
acc1 = mpu_raw(44:end, 4)';
acc2 = mtb_raw(30:end, 4)';
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

% get input model
[a_in, proc_var] = getModel(dt, horizon);
shift = 0.94;
a_in(1:length(a_in((shift/dt):end))) = a_in((shift/dt):end); % shift

% plot raw data
figure(1);
%plot(t, acc1, t, acc2, t, a_in);
plot(t, acc1, t, acc2);
legend("MPU Acc.", "MATLAB Acc.");
xlabel("Time (s)");
ylabel("Acceleration (m/s^2) - Gravity");
xlim([0 35]);

% integrate data
single_int = zeros(2, length(t));
z = zeros(2, length(t));
for i = 2:1:steps
    single_int(1, i) = single_int(1, i-1) + acc1(i-1)*dt;
    single_int(2, i) = single_int(2, i-1) + acc2(i-1)*dt;
    z(1, i) = z(1, i-1) + single_int(1, i-1)*dt;
    z(2, i) = z(2, i-1) + single_int(2, i-1)*dt;
end

% characterized noises
mpu_var = 4.4e-4;
mtb_var = 2.69e-5;
mtb_gm_var = 5.5e-4;
mtb_gm_beta = 0.3;
unity_wn = 1;
R = (1e-9)*eye(2); % set close to 0
% note: W must be calculated in loop due to changing proc_var

% system variables
n = 7;
A = [0 1 0 0 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 1 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 0 0 1 0;
     0 0 0 0 0 0 1;
     0 0 0 0 0 0 -mtb_gm_beta];
B = [0 0 0 0;
     1 0 0 0;
     0 0 0 0;
     0 1 0 0;
     0 0 0 0;
     0 0 1 0;
     0 0 0 sqrt(2*mtb_gm_var*mtb_gm_beta)];
G = [0; 1; 0; 0; 0; 0; 0];
H = [1 0 1 0 0 0 0;
     1 0 0 0 1 0 0];

% values of interest
default_pos = zeros(1, steps);
default_cov = zeros(1, steps);
obs_fix_pos = zeros(1, steps);
obs_fix_cov = zeros(1, steps);
low_proc_pos = zeros(1, steps);
low_proc_cov = zeros(1, steps);
low_proc_k = zeros(1, steps);
high_proc_pos = zeros(1, steps);
high_proc_cov = zeros(1, steps);
high_proc_k = zeros(1, steps);

% trial loop
% 1 -> Default procVar
% 2 -> Observability Fix
% 3 -> Low procVar
% 4 -> high proc var
for trial = 1:1:4

    % initial guesses
    xhatminus = zeros(n, 1);
    Pminus = zeros(n);
    
    % estimate and error placeholders
    xhat = zeros(n, steps);
    P = zeros(n, n, steps);

    % other placeholders
    sumK = zeros(1, steps);
    
    % Kalman filter loop
    for i = 1:1:steps
        
        % discretize with update proc var
        if (trial == 1 || trial == 2)
            W = eye(4) .* [proc_var(i), mpu_var, mtb_var, unity_wn];
        end
        if (trial == 3)
            W = eye(4) .* [1e-9, mpu_var, mtb_var, unity_wn];
        end
        if (trial == 4)
            W = eye(4) .* [1, mpu_var, mtb_var, unity_wn];
        end
        % van loan discretization
        L = [-A, B*W*B'; zeros(n), A'] * dt;
        E = expm(L);
        Phi = E(n+1:end, n+1:end)';
        Q = Phi * E(1:n, n+1:end);
    
        % step 1: compute Kalman gain
        K = Pminus * H' / (H*Pminus*H'+ R);
    
        % step 2: get measurement and update estimate
        xhat(:, i) = xhatminus + K*(z(:, i) - H*xhatminus);
    
        % step 3: update error covariance
        P(:, :, i) = (eye(n) - K*H) * Pminus;
    
        % step 4: project ahead
        % note the additional bit for acceleration input
        xhatminus = (Phi * xhat(:, i)) + (G*a_in(i)*dt);
        Pminus = Phi*P(:, :, i)*Phi' + Q;
    
        % use altenative model to "observe" states and inflection points.
        if (trial == 2)
            if ((t(i) >= 11.1) && (t(i) < 23.3)) || (t(i) >= 32.9)
                A_alt = [1 dt; 0 1];
                B_alt = [0; 0];
                C_alt = [1 0];
                D_alt = [1];
                On = [C_alt; C_alt*A_alt];
                Tn = [D_alt, zeros(1);
                      C_alt*B_alt, D_alt];
                true_bias1 = pinv(On) * ([z(1, i-1); z(1, i)] - (Tn*[xhatminus(1); xhatminus(1)]));
                true_bias2 = pinv(On) * ([z(2, i-1); z(2, i)]  - (Tn*[xhatminus(1); xhatminus(1)]));
                xhatminus(1) = xhatminus(1);
                xhatminus(2) = 0;
                xhatminus(3) = true_bias1(1); 
                xhatminus(4) = true_bias1(2); 
                xhatminus(5) = true_bias2(1); 
                xhatminus(6) = true_bias2(2); 
                Pminus = eye(7) .* [0, 0, 0, 0, 0, 0, Pminus(7,7)];
            end
        end
    
        % gather data
        sumK(i) = abs(K(1,1)) + abs(K(2,2));
    
    end
    
    if (trial == 1)
        default_pos = xhat(1, :);
        default_cov(:) = P(1, 1, :);
    elseif (trial == 2)
        obs_fix_pos = xhat(1, :);
        obs_fix_cov(:) = P(1, 1, :);
    elseif (trial == 3)
        low_proc_pos = xhat(1, :);
        low_proc_cov(:) = P(1, 1, :);
        low_proc_k = sumK;
    elseif (trial == 4)
        high_proc_pos = xhat(1, :);
        high_proc_cov(:) = P(1, 1, :);
        high_proc_k = sumK;
    end

end

% plots position and error covariance
for trial = 1:1:4
    figure(trial+1);
    subplot(211);
    if (trial == 1)
        plot(t, z(1, :), t, z(2, :), t, default_pos);
    elseif (trial == 2)
        plot(t, z(1, :), t, z(2, :), t, obs_fix_pos);
    elseif (trial == 3)
        plot(t, z(1, :), t, z(2, :), t, low_proc_pos);
    else
        plot(t, z(1, :), t, z(2, :), t, high_proc_pos);
    end
    xlabel("Time (s)"); ylabel('Position (m)');
    legend('MPU', 'MATLAB', 'Estimate');
    xlim([0 35]);
    subplot(212);
    plot(t, default_cov);
    if (trial == 1)
        plot(t, default_cov);
    elseif (trial == 2)
        plot(t, obs_fix_cov);
    elseif (trial == 3)
        plot(t, low_proc_cov);
    else
        plot(t, high_proc_cov);
    end
    xlabel("Time (s)"); ylabel('MSE');
    xlim([0 35]);
end

% compares kalman gains
figure(6);
plot(t, low_proc_k, t, high_proc_k);
xlabel("Time (s)"); ylabel('|K(1,1)| + |K(2,2)|');
legend('Low Proc. Var.', 'High Proc. Var');
xlim([-2 40]); ylim([-2 55]);

% plots all position estimators
figure(7);
plot(t, z(1, :), t, z(2, :), t, default_pos, t, obs_fix_pos, t, low_proc_pos, t, high_proc_pos);
xlabel("Time (s)"); ylabel('Position (m)');
xlim([0 35]);
legend('MPU', 'MATLAB', 'Default Estimate', 'Obs. Fix Estimate', 'Low Process Var. Estimate', 'High Process Var. Estimate');

% plots all covariances
figure(8);
plot(t, default_cov, t, obs_fix_cov, t, low_proc_cov, t, high_proc_cov);
xlabel("Time (s)"); ylabel('MSE');
xlim([0 35]);
legend('Default Estimate', 'Obs. Fix Estimate', 'Low Process Var. Estimate', 'High Process Var. Estimate');

% display heights and drifts
t1 = 10.5;
fprintf("%.2f & %.2f & %.2f & %.2f & %.2f & %.2f &  \\\\ \\hline \n", z(1, t1/dt), z(2, t1/dt), default_pos(t1/dt), obs_fix_pos(t1/dt), low_proc_pos(t1/dt), high_proc_pos(t1/dt));
t1 = 31.5;
fprintf("%.2f & %.2f & %.2f & %.2f & %.2f & %.2f &  \\\\ \\hline \n", z(1, t1/dt), z(2, t1/dt), default_pos(t1/dt), obs_fix_pos(t1/dt), low_proc_pos(t1/dt), high_proc_pos(t1/dt));