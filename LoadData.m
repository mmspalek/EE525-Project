clear;
clc;
close all;

tic
addpath('C:\Users\Precision\Desktop\EE525\Projects\Project2\Motion Data');



% Load data from the accelerometers
mat_test1 = load('C:\Users\Precision\Downloads\test1.mat');
mat_test2 = load('C:\Users\Precision\Downloads\test2.mat');
mat_test3 = load('C:\Users\Precision\Downloads\test3.mat');
mat_test4 = load('C:\Users\Precision\Downloads\test4.mat');
mat_test5 = load('C:\Users\Precision\Downloads\test5.mat');
mat_test6 = load('C:\Users\Precision\Downloads\test6.mat');
mat_test7 = load('C:\Users\Precision\Downloads\test7.mat');
mat_test8 = load('C:\Users\Precision\Downloads\test8.mat');
mat_test9 = load('C:\Users\Precision\Downloads\test9.mat');
mat_test10 = load('C:\Users\Precision\Downloads\test10.mat');
mat_test11 = load('C:\Users\Precision\Downloads\test11.mat');
mat_test12 = load('C:\Users\Precision\Downloads\test12.mat');
mat_test13 = load('C:\Users\Precision\Downloads\test13.mat');
mat_test14 = load('C:\Users\Precision\Downloads\test14.mat');
mat_test15 = load('C:\Users\Precision\Downloads\test15.mat');
mat_test16 = load('C:\Users\Precision\Downloads\test16.mat');


% Extracting the intended data from the accelerometer data
% Acceleration data from matlab sensor
test1 = mat_test1.Acceleration.Z;
test2 = mat_test2.Acceleration.Z;
test3 = mat_test3.Acceleration.Z;
test4 = mat_test4.Acceleration.Z;
test5 = mat_test5.Acceleration.Z;
test6 = mat_test6.Acceleration.Z;
test7 = mat_test7.Acceleration.Z;
test8 = mat_test8.Acceleration.Z;
test9 = mat_test9.Acceleration.Z;
test10 = mat_test10.Acceleration.Z;
test11 = mat_test11.Acceleration.Z;
test12 = mat_test12.Acceleration.Z;
test13 = mat_test13.Acceleration.Z;
test14 = mat_test14.Acceleration.Z;
test15 = mat_test15.Acceleration.Z;
test16 = mat_test16.Acceleration.Z;

%MPU Data
mpu_table = readtable('C:\Users\Precision\Desktop\EE525\Projects\Project3\Matlab\Copy of MPU6050_upDn.xlsx');
mpu_t4 = mpu_table{:,19};
mpu_t8 = mpu_table{:,23};
mpu_t9 = mpu_table{:,24};
mpu_t10 = mpu_table{:,25};
mpu_t15 = mpu_table{:,30};
% Selected 5 trial runs with similar elevator acceleration start time
% Stage 1 - Before accel
% Stage 2 - Travel to 2nd floor
% Stage 3 - Rest at 2nd floor
% Stage 4 - Travel down to 1st
% Stage 5 - Rest at 1st

f = 50;
% % Extraction of 5 different stages from matlab and mpu data
% test15_s1 = test15(1:3.28*f+1);
% test15_s2 = test15(3.28*f+1:12.1*f+1);
% test15_s3 = test15(12.1*f+1:24.3*f+1);
% test15_s4 = test15(24.3*f+1:33.9*f+1);
% test15_s5 = test15(33.9*f+1:41*f+1);
% 
% test14_s1 = test14(1:3.28*f+1);
% test14_s2 = test14(3.28*f+1:12.2*f+1);
% test14_s3 = test14(12.2*f+1:24.2*f+1);
% test14_s4 = test14(24.2*f+1:33.7*f+1);
% test14_s5 = test14(33.7*f+1:41*f+1);
% 
% test13_s1 = test13(1:3.34*f+1);
% test13_s2 = test13(3.34*f+1:12.4*f+1);
% test13_s3 = test13(12.4*f+1:24.5*f+1);
% test13_s4 = test13(24.5*f+1:34*f+1);
% test13_s5 = test13(34*f+1:41*f+1);
% 
% test8_s1 = test8(1:3.26*f+1);
% test8_s2 = test8(3.26*f+1:12.3*f+1);
% test8_s3 = test8(12.3*f+1:24.7*f+1);
% test8_s4 = test8(24.7*f+1:34.4*f+1);
% test8_s5 = test8(34.4*f+1:41*f+1);
% 
% test5_s1 = test5(1:3.2*f+1);
% test5_s2 = test5(3.2*f+1:12.4*f+1);
% test5_s3 = test5(12.4*f+1:24.6*f+1);
% test5_s4 = test5(24.6*f+1:34.5*f+1);
% test5_s5 = test5(34.5*f+1:41*f+1);
% 
% mput4_s1 = mpu_t4(1:2.98*f+1);
% mput4_s2 = mpu_t4(2.98*f+1:12.3*f+1);
% mput4_s3 = mpu_t4(12.3*f+1:24.56*f+1);
% mput4_s4 = mpu_t4(24.56*f+1:34*f+1);
% mput4_s5 = mpu_t4(34*f+1:41*f+1);
% 
% mput8_s1 = mpu_t8(1:3.1*f+1);
% mput8_s2 = mpu_t8(3.1*f+1:12.3*f+1);
% mput8_s3 = mpu_t8(12.3*f+1:24.5*f+1);
% mput8_s4 = mpu_t8(24.5*f+1:33.7*f+1);
% mput8_s5 = mpu_t8(33.7*f+1:41*f+1);
% 
% mput9_s1 = mpu_t9(1:2.92*f+1);
% mput9_s2 = mpu_t9(2.92*f+1:12.1*f+1);
% mput9_s3 = mpu_t9(12.1*f+1:23.9*f+1);
% mput9_s4 = mpu_t9(23.9*f+1:33.6*f+1);
% mput9_s5 = mpu_t9(33.6*f+1:41*f+1);
% 
% mput10_s1 = mpu_t4(1:3*f+1);
% mput10_s2 = mpu_t4(3*f+1:12.2*f+1);
% mput10_s3 = mpu_t4(12.2*f+1:24.2*f+1);
% mput10_s4 = mpu_t4(24.2*f+1:34*f+1);
% mput10_s5 = mpu_t4(34*f+1:41*f+1);
% 
% mput15_s1 = mpu_t4(1:3*f+1);
% mput15_s2 = mpu_t4(3*f+1:11.8*f+1);
% mput15_s3 = mpu_t4(11.8*f+1:24*f+1);
% mput15_s4 = mpu_t4(24*f+1:33.5*f+1);
% mput15_s5 = mpu_t4(33.5*f+1:41*f+1);
% 
% % Variance calculations
% var15_s1 = var(test15_s1); var15_s2 = var(test15_s2); var15_s3 = var(test15_s3); var15_s4 = var(test15_s4); var15_s5 = var(test15_s5);
% 
% var14_s1 = var(test14_s1); var14_s2 = var(test14_s2); var14_s3 = var(test14_s3); var14_s4 = var(test14_s4); var14_s5 = var(test14_s5);
% 
% var13_s1 = var(test13_s1); var13_s2 = var(test13_s2); var13_s3 = var(test13_s3); var13_s4 = var(test13_s4); var13_s5 = var(test13_s5);
% 
% var8_s1 = var(test8_s1); var8_s2 = var(test8_s2); var8_s3 = var(test8_s3); var8_s4 = var(test8_s4); var8_s5 = var(test8_s5);
% 
% var5_s1 = var(test5_s1); var5_s2 = var(test5_s2); var5_s3 = var(test5_s3); var5_s4 = var(test5_s4); var5_s5 = var(test5_s5);
% 
% 
% mpuvar15_s1 = var(mput15_s1); mpuvar15_s2 = var(mput15_s2); mpuvar15_s3 = var(mput15_s3); mpuvar15_s4 = var(mput15_s4); mpuvar15_s5 = var(mput15_s5);
% 
% mpuvar10_s1 = var(mput10_s1); mpuvar10_s2 = var(mput10_s2); mpuvar10_s3 = var(mput10_s3); mpuvar10_s4 = var(mput10_s4); mpuvar10_s5 = var(mput10_s5);
% 
% mpuvar9_s1 = var(mput9_s1); mpuvar9_s2 = var(mput9_s2); mpuvar9_s3 = var(mput9_s3); mpuvar9_s4 = var(mput9_s4); mpuvar9_s5 = var(mput9_s5);
% 
% mpuvar8_s1 = var(mput8_s1); mpuvar8_s2 = var(mput8_s2); mpuvar8_s3 = var(mput8_s3); mpuvar8_s4 = var(mput8_s4); mpuvar8_s5 = var(mput8_s5);
% 
% mpuvar4_s1 = var(mput4_s1); mpuvar4_s2 = var(mput4_s2); mpuvar4_s3 = var(mput4_s3); mpuvar4_s4 = var(mput4_s4); mpuvar4_s5 = var(mput4_s5);
% 
% % Variance values for each dataset
% var15 = [var15_s1, var15_s2, var15_s3, var15_s4, var15_s5, mpuvar15_s1, mpuvar15_s2, mpuvar15_s3, mpuvar15_s4, mpuvar15_s5];
% var14 = [var14_s1, var14_s2, var14_s3, var14_s4, var14_s5, mpuvar10_s1, mpuvar10_s2, mpuvar10_s3, mpuvar10_s4, mpuvar10_s5];
% var13 = [var13_s1, var13_s2, var13_s3, var13_s4, var13_s5, mpuvar9_s1, mpuvar9_s2, mpuvar9_s3, mpuvar9_s4, mpuvar9_s5];
% var8  = [var8_s1, var8_s2, var8_s3, var8_s4, var8_s5, mpuvar8_s1, mpuvar8_s2, mpuvar8_s3, mpuvar8_s4, mpuvar8_s5];
% var5  = [var5_s1, var5_s2, var5_s3, var5_s4, var5_s5, mpuvar4_s1, mpuvar4_s2, mpuvar4_s3, mpuvar4_s4, mpuvar4_s5];
% 
% % Combine all data into a single matrix (rows: trials, columns: stages)
% Data = [var15; var14; var13; var8; var5];
% 
% % Create a new table with switched rows and columns
% VarianceTable = array2table(Data', ...
%     'VariableNames', {'Variance Trial1', 'Variance Trial2', 'Variance Trial3', 'Variance Trial4', 'Variance Trial5'}, ...
%     'RowNames', {'Stage 1 MAT', 'Stage 2 MAT', 'Stage 3 MAT', 'Stage 4 MAT', 'Stage 5 MAT', 'Stage 1 MPU', 'Stage 2 MPU', 'Stage 3 MPU', 'Stage 4 MPU', 'Stage 5 MPU'});
% 
% % Display the transposed table
% disp(VarianceTable);

%% More data plotting
% system variables
T = 0.02;

% read-in data and get average
mpu_table = readtable('MPU6050_upDn.xlsx');
mtb_table = readtable('MATLAB_upDn.xlsx');
mpu_raw = rmmissing(mpu_table{:, 1:15}) * (9.81 / 2^14);
mat_raw = rmmissing(mtb_table{:, 1:15});
avg = (mean(mpu_raw, 2) + mean(mat_raw, 2)) * 0.5;

mpu_cln = zeros(size(mpu_raw));
mat_cln = zeros(size(mat_raw));

mpu_cln(1:1000,:) = mpu_raw(1:1000,:) - mean(mpu_raw(1:2/T,:),1); % remove mean of first 2 seconds
mpu_cln(1001:end,:) = mpu_raw(1001:end,:) - mean(mpu_raw((34.5/T):(36.5/T),:),1); % remove mean of first 2 seconds
mat_cln(1:1000,:) = mat_raw(1:1000,:) - mean(mat_raw(1:2/T,:) + 0.001,1); % remove mean of first 2 seconds
mat_cln(1001:end,:) = mat_raw(1001:end,:) - mean(mat_raw((34.5/T):(36.5/T),:) +0.001,1); % remove mean of first 2 seconds
tic

% mpu_cln = mpu_raw - mean(mpu_raw(1:2/T));

t = (0:T:length(mat_cln)*T-T);
% Velocity & Height for MPU
for j = 1:width(mpu_cln)
mpu_v(:,j) = cumtrapz(t,mpu_cln(:,j));
mpu_p(:,j) = cumtrapz(t,mpu_v(:,j));
end

% Velocity & Height for MAT
for j = 1:width(mat_cln)
mat_v(:,j) = cumtrapz(t,mat_cln(:,j));
mat_p(:,j) = cumtrapz(t,mat_v(:,j));
end

tbuff(:,1) = [22; 81; 62; 44; 56; 82; 65; 49; 39; 53; 37; 77; 69; 72; 43];
tbuff(:,2) = [1; 53; 44; 30; 59; 41; 48; 70; 54; 42; 48; 29; 61; 60; 58];


%% Plot of Height for MPU Trials
colors = lines(15); % Use 15 distinct colors

figure;
hold on;

for i = 1:size(tbuff, 1)
    plot(mpu_p(tbuff(i, 1):end, i), 'Color', colors(i, :), 'LineStyle', '-');
end
yline(0, '--', 'linewidth', 0.5)
title('First 15 MPU trials')
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'Trial 6', 'Trial 7', 'Trial 8', 'Trial 9', 'Trial 10', 'Trial 11', 'Trial 12', 'Trial 13', 'Trial 14', 'Trial 15')

% Plot of Height for MAT Trials
figure;
hold on;
for i = 1:size(tbuff, 1)
    plot(mat_p(tbuff(i, 2):end, i), 'Color', colors(i, :), 'LineStyle', '-');
end
yline(0, '--', 'linewidth', 0.5)
title('First 15 MAT trials')
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'Trial 6', 'Trial 7', 'Trial 8', 'Trial 9', 'Trial 10', 'Trial 11', 'Trial 12', 'Trial 13', 'Trial 14', 'Trial 15')

% hcolumn_mpu = mpu_p(end,:)
% hcolumn_mat = mat_p(end,:)

%%

% Determing height at different inflection pts
column_vector_mpu = []; % Initialize an empty column vector for mpu_p
column_vector_mat = []; % Initialize an empty column vector for mat_p
disp('Height at 1st');
for j = 1:width(mpu_cln)
    column_vector_mpu = [column_vector_mpu; mpu_p(tbuff(j,1), j)];
    column_vector_mat = [column_vector_mat; mat_p(tbuff(j,2), j)];
end
% Display the column vectors
disp('Column vector at 1st for mpu_p:');
disp(column_vector_mpu);
disp('Column vector at 1st for mat_p:');
disp(column_vector_mat);


% Height at 1st inflection
column_vector_mpu = []; % Initialize an empty column vector for mpu_p
column_vector_mat = []; % Initialize an empty column vector for mat_p
disp('Height arriving at 1st');
for j = 1:width(mpu_cln)
    column_vector_mpu = [column_vector_mpu; mpu_p(500+tbuff(j,1), j)];
    column_vector_mat = [column_vector_mat; mat_p(500+tbuff(j,2), j)];
end
% Display the column vectors
disp('Column vector at 2nd for mpu_p:');
disp(column_vector_mpu);
disp('Column vector at 2nd for mat_p:');
disp(column_vector_mat);


% Height at 2nd inflection pt
column_vector_mpu = []; % Initialize an empty column vector for mpu_p
column_vector_mat = []; % Initialize an empty column vector for mat_p
disp('Height arriving leaving 2nd');
for j = 1:width(mpu_cln)
    column_vector_mpu = [column_vector_mpu; mpu_p(1200+tbuff(j,1), j)];
    column_vector_mat = [column_vector_mat; mat_p(1200+tbuff(j,2), j)];
end
% Display the column vectors
disp('Column vector leaving 2nd for mpu_p:');
disp(column_vector_mpu);
disp('Column vector leaving 2nd for mat_p:');
disp(column_vector_mat);


% Height at 3rd inflection pt
column_vector_mpu = []; % Initialize an empty column vector for mpu_p
column_vector_mat = []; % Initialize an empty column vector for mat_p

disp('Height arriving leaving 2nd');
for j = 1:width(mpu_cln)
    column_vector_mpu = [column_vector_mpu; mpu_p(1550+tbuff(j,1), j)];
    column_vector_mat = [column_vector_mat; mat_p(1550+tbuff(j,2), j)];
end
% Display the column vectors
disp('Column vector back at 1st for mpu_p:');
disp(column_vector_mpu);
disp('Column vector back at 1st for mat_p:');
disp(column_vector_mat);

%%
% Plot of OG Data from MPU
figure;
hold on;

for i = 1:size(tbuff, 1)
    t = (tbuff(i,1):length(mpu_cln(:,i)))/50;
    plot(t,mpu_cln(tbuff(i, 1):end, i));
end


title('First 15 MPU trials')
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'Trial 6', 'Trial 7', 'Trial 8', 'Trial 9', 'Trial 10', 'Trial 11', 'Trial 12', 'Trial 13', 'Trial 14', 'Trial 15')
ylim([-1 1])
xlim([0 40])
ylabel('Acceleration (m/s^2)')
xlabel('Time (s)')
grid on; box;
% Plot of OG Data from MAT
figure;
hold on;

for i = 1:size(tbuff,1)
     t = (tbuff(i,2):length(mpu_cln(:,i)))/50;
    plot(t,mat_cln(tbuff(i,2):end, i));
end

title('First 15 MATLAB trials')
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'Trial 6', 'Trial 7', 'Trial 8', 'Trial 9', 'Trial 10', 'Trial 11', 'Trial 12', 'Trial 13', 'Trial 14', 'Trial 15')
ylabel('Acceleration (m/s^2)')
xlabel('Time (s)')
grid on; box;
ylim([-1 1])

%% Plot of best trials
matt = [4, 8, 11];
mput = [2, 3, 9];

tbuff(2,1) = 30;
% tbuff(4,2) = tbuff(4,2) + 20;

figure;
subplot(211)

hold on;
t = (tbuff(mput(1),1)-30:length(mpu_cln(:,mput(1)))-30)/50;
plot(t, mpu_cln(tbuff(mput(1),1):end,mput(1)), 'DisplayName', sprintf('Trial %d', 1));
t = (tbuff(mput(2),1):length(mpu_cln(:,mput(2))))/50;
plot(t, mpu_cln(tbuff(mput(2),1):end,mput(2)), 'DisplayName', sprintf('Trial %d', 2));
t = (tbuff(mput(3),1):length(mpu_cln(:,mput(3))))/50;
plot(t, mpu_cln(tbuff(mput(3),1):end,mput(3)), 'DisplayName', sprintf('Trial %d', 3));

title('Acceleration data for different MPU Trials');
xlabel('Time');
ylabel('Acceleration (m/s^2)');
legend show; % Automatically display legend with trial labels
grid on;
ylim([-0.7 0.7])

subplot(212)

hold on;
i = 1;
t = (tbuff(matt(i),2)+36:length(mat_cln(:,matt(i)))+36)/50;
plot(t, mat_cln(tbuff(matt(i),2):end,matt(i)), 'DisplayName', sprintf('Trial %d', 1));
i = 2;
t = (tbuff(matt(i),2):length(mat_cln(:,matt(i))))/50;
plot(t, mat_cln(tbuff(matt(i),2):end,matt(i)), 'DisplayName', sprintf('Trial %d', 2));
i = 3;
t = (tbuff(matt(i),2):length(mat_cln(:,matt(i))))/50;
plot(t, mat_cln(tbuff(matt(i),2):end,matt(i)), 'DisplayName', sprintf('Trial %d', 3));

title('Acceleration data for different MAT Trials');
xlabel('Time');
ylabel('Acceleration (m/s^2)');
legend show; % Automatically display legend with trial labels
grid on;
ylim([-0.7 0.7])


%%
% Velocity Data
figure;
subplot(211)

hold on;
t = (tbuff(mput(1),1)-30:length(mpu_cln(:,mput(1)))-30)/50;
plot(t, mpu_v(tbuff(mput(1),1):end,mput(1)), 'DisplayName', sprintf('Trial %d', 1), 'linewidth', 1.5);
t = (tbuff(mput(2),1):length(mpu_cln(:,mput(2))))/50;
plot(t, mpu_v(tbuff(mput(2),1):end,mput(2)), 'DisplayName', sprintf('Trial %d', 2), 'linewidth', 1.5);
t = (tbuff(mput(3),1):length(mpu_cln(:,mput(3))))/50;
plot(t, mpu_v(tbuff(mput(3),1):end,mput(3)), 'DisplayName', sprintf('Trial %d', 3), 'linewidth', 1.5);

title('Velocity data for different MPU Trials');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend show; % Automatically display legend with trial labels
grid on; box
ylim([-0.6 0.6])
xlim([0 40])


subplot(212)
hold on;
i = 1;
t = (tbuff(matt(i),2)+36:length(mat_cln(:,matt(i)))+36)/50;
plot(t, mat_v(tbuff(matt(i),2):end,matt(i)), 'DisplayName', sprintf('Trial %d', 1), 'linewidth', 1.5);
i = 2;
t = (tbuff(matt(i),2):length(mat_cln(:,matt(i))))/50;
plot(t, mat_v(tbuff(matt(i),2):end,matt(i)), 'DisplayName', sprintf('Trial %d', 2), 'linewidth', 1.5);
i = 3;
t = (tbuff(matt(i),2):length(mat_cln(:,matt(i))))/50;
plot(t, mat_v(tbuff(matt(i),2):end,matt(i)), 'DisplayName', sprintf('Trial %d', 3), 'linewidth', 1.5);

title('Velocity data for different MAT Trials');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend show; % Automatically display legend with trial labels
grid on; box
ylim([-0.6 0.6])
xlim([0 40])

%%
% Height Data
figure;

% Subplot 1: MPU Trials
subplot(211)
hold on;

% Plot MPU Trials
t = (tbuff(mput(1),1)-30:length(mpu_cln(:,mput(1)))-30)/50;
plot(t, mpu_p(tbuff(mput(1),1):end,mput(1)), 'DisplayName', sprintf('Trial %d', 1), 'linewidth', 1.5);
t = (tbuff(mput(2),1):length(mpu_cln(:,mput(2))))/50;
plot(t, mpu_p(tbuff(mput(2),1):end,mput(2)), 'DisplayName', sprintf('Trial %d', 2), 'linewidth', 1.5);
t = (tbuff(mput(3),1):length(mpu_cln(:,mput(3))))/50;
plot(t, mpu_p(tbuff(mput(3),1):end,mput(3)), 'DisplayName', sprintf('Trial %d', 3), 'linewidth', 1.5);

% Title, Labels, and Formatting
title('Height data for different MPU Trials');
xlabel('Time (s)');
ylabel('Height (m)');
legend show; % Automatically display legend with trial labels
grid on; box

% Add yline with HandleVisibility set to 'off'
yline(0, 'HandleVisibility', 'off');
ylim([-3 5]);
xlim([0 40]);

% Subplot 2: MAT Trials
subplot(212)
hold on;

% Plot MAT Trials
i = 1;
t = (tbuff(matt(i),2)+36:length(mat_cln(:,matt(i)))+36)/50;
plot(t, mat_p(tbuff(matt(i),2):end,matt(i)), 'DisplayName', sprintf('Trial %d', 1), 'linewidth', 1.5);
i = 2;
t = (tbuff(matt(i),2):length(mat_cln(:,matt(i))))/50;
plot(t, mat_p(tbuff(matt(i),2):end,matt(i)), 'DisplayName', sprintf('Trial %d', 2), 'linewidth', 1.5);
i = 3;
t = (tbuff(matt(i),2):length(mat_cln(:,matt(i))))/50;
plot(t, mat_p(tbuff(matt(i),2):end,matt(i)), 'DisplayName', sprintf('Trial %d', 3), 'linewidth', 1.5);

% Title, Labels, and Formatting
title('Height data for different MAT Trials');
xlabel('Time (s)');
ylabel('Height (m)');
legend show; % Automatically display legend with trial labels
grid on; box

% Add yline with HandleVisibility set to 'off'
yline(0, 'HandleVisibility', 'off');
ylim([-3 5]);
xlim([0 40]);


%%
% Calculate variances for each segment and store them as variables
for i = 1:size(tbuff, 1)
    var_s1(i) = var(mpu_cln(tbuff(i, 1):tbuff(i, 1) + 110, i));
    var_s2(i) = var(mpu_cln(tbuff(i, 1) + 110:tbuff(i, 1) + 250, i));
    var_s3(i) = var(mpu_cln(tbuff(i, 1) + 250:tbuff(i, 1) + 450, i));
    var_s4(i) = var(mpu_cln(tbuff(i, 1) + 450:tbuff(i, 1) + 525, i));
    var_s5(i) = var(mpu_cln(tbuff(i, 1) + 525:tbuff(i, 1) + 1160, i));
    var_s6(i) = var(mpu_cln(tbuff(i, 1) + 1160:tbuff(i, 1) + 1260, i));
end

% Combine the variables into a matrix
var_matrix = [var_s1(:), var_s2(:), var_s3(:), var_s4(:), var_s5(:), var_s6(:)];

% Compute the column averages
avg_row = mean(var_matrix, 1);

% Create a table with rows for trials and columns for segments
var_table = array2table(var_matrix, ...
    'VariableNames', {'Stage1', 'Stage2', 'Stage3', 'Stage4', 'Stage5', 'Stage6'}, ...
    'RowNames', strcat("Trial_", string(1:size(tbuff, 1))));

% Add the averages as the final row
var_table{end + 1, :} = avg_row;
var_table.Properties.RowNames{end} = 'Average';

% Display the table
disp(var_table);


%% State space model(s)
A = [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
B = [0 0 0; 0 1 0; 0 0 0; 0 0 1];
C = [1 0 0 0; 0 0 1 0];
D = zeros(2,3);
sys = ss(A,B,C,D);



observability = rank(obsv(A, C))

if observability == height(A)
    display('This model is fully observable');
else 
    display('This model is not fully observable');
end


[V, lam] = eig(A)


%%
% Last value variance(s)
mpulastval(1) = mpu_p(end,1);
mpulastval(2) = mpu_p(end,2);
mpulastval(3) = mpu_p(end,3);
mpulastval(4) = mpu_p(end,9);
mpulastval(5) = mpu_p(end,11);

lastval_mean_mpu = mean(mpulastval);
lastval_var_mpu = var(mpulastval);

matlastval(1) = mat_p(end,4);
matlastval(2) = mat_p(end,6);
matlastval(3) = mat_p(end,8);
matlastval(4) = mat_p(end,5);
matlastval(5) = mat_p(end,12);

lastval_mean_mat = mean(matlastval)
lastval_var_mat = var(matlastval)