% Nolan Egging
% EE525 - Draft 3
% Data Analysis
% Due Nov 22

close all;
clear;
clc;

% system variables
T = 1;

% state space model 1
A_ct = [0 1 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 1 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 1;
        0 0 0 0 0 0];
A = expm(A_ct*T);
C = [1 0 1 0 0 0;
     1 0 0 0 1 0];
On = [C; C*A; C*A*A];
disp("Rank of SS 1: " + rank(On));
G = On'*On;
[V, L] = eig(G);

%%

% state space model 2
A_ct = [0 1;
        0 0];
A = expm(A_ct*T);
C = [1 0;
     0 0];
On = [C; C*A];
disp("Rank of SS 2: " + rank(On));
G = On'*On;
[V, L] = eig(G);

% plot observability ellipsoids
v1 = L(1,1)^(-0.5) * V(:, 1);
v2 = L(2,2)^(-0.5) * V(:, 2);
[t1, t2] = createEllipse(v1, v2);
figure(2);
plot(t1, t2, [0, v1(1)], [0, v1(2)], [0, v2(1)], [0, v2(2)]);
xlabel("Bias");
ylabel("Derivative of Bias");
legend("Observability Ellispoid", "Eigenvector 1 * 1/sqrt(\lambda)", "Eigenvector 2 * 1/sqrt(\lambda)")
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);

% equally observable
theta = -angle(v1(1) + sqrt(-1)*v1(2));
Tf = [L(1,1)^0.5, 0; 0, L(2,2)^0.5] * [cos(theta), -sin(theta); sin(theta), cos(theta)];
A = Tf*A*inv(Tf);
C = C*inv(Tf);
On = [C; C*A];
G = On'*On;
[V, L] = eig(G);

% plot equally observability ellipsoids
v1 = L(1,1)^(-0.5) * V(:, 1);
v2 = L(2,2)^(-0.5) * V(:, 2);
[t1, t2] = createEllipse(v1, v2);
figure(3);
plot(t1, t2, [0, v1(1)], [0, v1(2)], [0, v2(1)], [0, v2(2)]);
xlabel("Bias");
ylabel("Derivative of Bias");
legend("Observability Ellispoid", "Eigenvector 1 * 1/sqrt(\lambda)", "Eigenvector 2 * 1/sqrt(\lambda)")
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);

% creates ellipse given two orthogonal vectors.
function [x, y] = createEllipse(vec1, vec2)
    % assumes vectors are orthogonal

    j = sqrt(-1);
    theta = 0:pi/1800:2*pi;

    a = vec1(1) + j*vec1(2);
    b = vec2(1) + j*vec2(2);
    rotation = angle(a);

    pre_x = abs(a) * cos(theta);
    pre_y = abs(b) * sin(theta);
    x = pre_x*cos(rotation) - pre_y*sin(rotation);
    y = pre_y*cos(rotation) + pre_x*sin(rotation);
    
end