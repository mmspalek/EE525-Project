% Nolan Egging
% EE525 - Draft 3
% Model Verification
% Due Nov 22

close all;
clear;
clc;

% system variables
T = 0.02;
horizon = 40; 
t = 0:T:horizon;

% noise models
proc_var = 0;
mpu_var = 0.00043538;
matlab_wn_var = 0.0000269;
matlab_cor_var = 0.0000055;
matlab_beta = 0.3;

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

% create input vectors
u = zeros(3, length(t));
[a, p] = getModel(T, horizon);
u(1, :) = a;

% no process noise
w = zeros(height(A), length(t));

% output and state vectors with no noise
[x, y] = applyDTModel(A, B, C, D, u, T*w);

% plot inital results
figure(1);
subplot(221);
plot(t, x(1, :));
title("State Vectors: Displacement");
ylim([-0.2 3.5]);
subplot(222);
plot(t, x(3, :), t, x(5, :));
title("State Vectors: Bias");
legend("b1", "b2");
ylim([-0.1 0.1]);
subplot(223);
plot(t, x(4, :), t, x(6, :));
title("State Vectors: Derivative of Bias");
legend("b1", "b2");
ylim([-0.1 0.1]);
subplot(224);
plot(t, y(1, :), t, y(2, :));
title("Output Vectors")
lg1 = legend("y1", "y2");
lg1.Location = "south";
ylim([-0.2 3.5]);

% create random noise inputs
rng(1);
u(2, :) = sqrt(mpu_var)*randn(1, length(t));
u(3, :) = createGMNoise(matlab_cor_var, matlab_beta, T, t) + sqrt(matlab_wn_var)*randn(1, length(t));
[x, y] = applyDTModel(A, B, C, D, u, T*w);

% plot single with bias results
figure(2);
subplot(221);
plot(t, x(1, :));
title("State Vectors: Displacement");
ylim([-0.2 3.5]);
subplot(222);
plot(t, x(3, :), t, x(5, :));
title("State Vectors: Bias");
legend("b1", "b2");
subplot(223);
plot(t, x(4, :), t, x(6, :));
title("State Vectors: Derivative of Bias");
legend("b1", "b2");
subplot(224);
plot(t, y(1, :), t, y(2, :));
title("Output Vectors")
lg1 = legend("y1", "y2");
lg1.Location = "south";


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
G = Ok'*Ok;
[V, L] = eig(G);

% plot pseudo-inverse results
figure(3);
subplot(221);
hold on;
plot(t, x(1, :));
plot(t, xls(1, :));
hold off;
legend("Height", "LS Estimate");
subplot(222);
hold on;
plot(t, x(3, :));
plot(t, xls(3, :));
hold off;
legend("Bias 1", "LS Estimate");
subplot(223);
hold on;
plot(t, x(5, :));
plot(t, xls(5, :));
hold off;
legend("Bias 2", "LS Estimate");
subplot(224);
hold on;
plot(t, y(1, :));
plot(t, y(2, :));
plot(t, xls(1, :));
hold off;
legend("Meas. Height 1", "Meas. Height 2", "LS Estimate");

% run multiple sims to get averages and variances
sims = 100;
sum_x = zeros(6, length(t));
sum_y = zeros(2, length(t));
sum_xls = zeros(6, length(t));
final_y1 = zeros(1, sims);
final_y2 = zeros(1, sims);
final_xls = zeros(1, sims);
for j = 1:1:sims
    disp(j);

    u(2, :) = sqrt(mpu_var)*randn(1, length(t));
    u(3, :) = createGMNoise(matlab_cor_var, matlab_beta, T, t) + sqrt(matlab_wn_var)*randn(1, length(t));
    [x, y] = applyDTModel(A, B, C, D, u, T*w);

    % apply pseudo inverse
    for h = 1:1:length(t)
        if h+k-1 == length(t)
            break;
        end
        uk = [u(1, h:h+k-1); zeros(1, k); zeros(1, k)];
        uk = uk(:);
        yk = y(:, h:h+k-1);
        yk = yk(:);
        xls(:, h) = pinv(Ok) * (yk - (Tk*uk));
    end

    % add to sum for average
    sum_x = sum_x + x;
    sum_y = sum_y + y;
    sum_xls = sum_xls + xls;

    % gets "final" y1 and y2 for variance
    final_y1(j) = y(1, 33/T);
    final_y2(j) = y(2, 33/T);
    final_xls(j) = xls(1, 33/T);

end

% get average
avg_x = sum_x / sims;
avg_y = sum_y / sims;
avg_xls = sum_xls / sims;

% plot average results with bias
figure(4);
subplot(221);
plot(t, avg_x(1, :));
title("State Vectors: Displacement");
ylim([-0.2 3.5]);
subplot(222);
plot(t, avg_x(3, :), t, avg_x(5, :));
title("State Vectors: Bias");
legend("b1", "b2");
ylim([-0.05 0.05]);
subplot(223);
plot(t, avg_x(4, :), t, avg_x(6, :));
title("State Vectors: Derivative of Bias");
legend("b1", "b2");
ylim([-0.05 0.05]);
subplot(224);
plot(t, avg_y(1, :), t, avg_y(2, :));
title("Output Vectors")
lg1 = legend("y1", "y2");
lg1.Location = "south";
ylim([-0.2 3.5]);

% plot average least-square estimate
figure(5);
plot(t, x(1, :));
hold on;
plot(t, avg_xls(1, :));
hold off;
legend("Height", "Average LS Estimate");

% display variances
disp("Variance in y1 (at inflection at 33s): " + var(final_y1));
disp("Variance in y2 (at inflection at 33s): " + var(final_y2));
disp("Estimate of Least Squares Final Value (at inflection at 33s): " + var(final_xls));

% functions
function [w] = createGMNoise(var, beta, T, t)
    walk = sqrt(var)*sqrt(1-exp(-2*beta*T))*randn(1, length(t));
    w = zeros(1, length(t));
    w(1) = sqrt(var) * randn(1);
    for i = 2:1:length(t)
        w(i) = w(i-1) * exp(-1*beta*T) + walk(i);
    end
end