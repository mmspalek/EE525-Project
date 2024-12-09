% Nolan Egging
% Max Spalek (edited)
% EE525 - Draft 3
% Input Models
% Due Nov 22

close all;
clear;
clc;
tic

% system variables
T = 0.02;

% read-in data and get average
mpu_table = readtable('MPU6050_upDn.xlsx');
mtb_table = readtable('MATLAB_upDn.xlsx');
mpu_raw = rmmissing(mpu_table{:, 1:15}) * (9.81 / 2^14);
mtb_raw = rmmissing(mtb_table{:, 1:15});
avg = (mean(mpu_raw, 2) + mean(mtb_raw, 2)) * 0.5;

t = (0:T:(length(mpu_raw)-1)*T);

% clean up data
rmavg_data = zeros(1, length(avg));
rmavg_data(1:1000) = avg(1:1000) - mean(avg(1:2/T)); % remove mean of first 2 seconds
rmavg_data(1001:end) = avg(1001:end) - mean(avg((34.5/T):(36.5/T))); % remove mean of first 2 seconds
clean_data = rmavg_data(1:1:length(t));

% create model
model = getModel(T, (length(clean_data)-1)*T);

% plot model against mpu
figure(2);
plot(t, clean_data);
hold on;
plot(t, model(1:length(t)));
% plot(t,model_int)

hold off;
lgd1 = legend("30 Run Average", "Model", "Model Vel");
lgd1.Location = "southeast";
title("Sensor Acceleration versus Model");
xlabel("Time (s)");
ylabel("Input Acceleration (m/s/s)");
xlim([0 40]);




% integrate model
model_int = (0:T:(length(clean_data)-1)*T);
model_double = (0:T:(length(clean_data)-1)*T);
for i = 2:1:length(model_int)
    model_int(i) = model_int(i-1) + model(i-1)*T;
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
%% -11/19 work
clean_val = (0:T:(length(clean_data)-1)*T);
clean_pos = (0:T:(length(clean_data)-1)*T);
for i = 2:1:length(clean_val)
    clean_val(i) = clean_val(i-1) + clean_data(i-1)*T;
    clean_pos(i) = clean_pos(i-1) + clean_val(i-1)*T;
end

figure(4);
% First subplot: Velocity with color-coded segments
subplot(211);
hold on; % Keep multiple segments on the same plot
for i = 1:length(t)-1
    if clean_val(i) >= 0
        plot(t(i:i+1), clean_val(i:i+1), 'g', 'LineWidth', 2); % Green for positive
    else
        plot(t(i:i+1), clean_val(i:i+1), 'r', 'LineWidth', 2); % Red for negative
    end
end
title('Integral of clean data: Data Velocity');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
ylim([-0.6 0.6]);
hold off;

% Second subplot: Regular position plot
subplot(212);
plot(t, clean_pos);
title('Double Integral of Data: Data Position');
xlabel('Time (s)');
ylabel('Position (m)');
xlim([0 40]);
ylim([-0.5 4.0]);


figure(5);
plot(t, clean_data - model(1:length(t)))
hold on;
xline(130*T);xline(300*T);xline(483*T);xline(605*T); xline(817*T); xline(950*T);
xline(1200*T); xline(1331*T);xline(1531*T);xline(1652*T); xline(1838*T); xline(1925*T);
title('Difference of Model & Average of 30 Trials')
xlabel('Time (s)')

%% Plot of Model subtracted from data
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


%Plot of the 30 diff trials

figure(1);
tbuff(:,1) = [22; 81; 62; 44; 56; 82; 65; 49; 39; 53; 37; 77; 69; 72; 43];
tbuff(:,2) = [1; 53; 44; 30; 59; 41; 48; 70; 54; 42; 48; 29; 61; 60; 58];
hold on;
for i = 1:size(tbuff, 1)
    plot(mpu_cln(tbuff(i, 1):end, i));
end

for i = 1:size(tbuff,2)
    plot(mat_cln(tbuff(i,2):end, i));
end
title('All 30 trials - both sensors')
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'Trial 6', 'Trial 7', 'Trial 8', 'Trial 9', 'Trial 10', 'Trial 11', 'Trial 12', 'Trial 13', 'Trial 14', 'Trial 15')



mpu_diff_cell = cell(1, 15);
mat_diff_cell = mpu_diff_cell;
 %Plot of Difference of Model & Data
% for i = 1:length(tbuff)
% 
%     figure;
%     subplot(211) % MPU Diff plots
%     mpu_diff =  mpu_cln(tbuff(i,1):end,i)' - model(1:length(mpu_cln(:,i)) - tbuff(i,1)+1);
%     mpu_diff_cell{i} = mpu_diff; % Store the vector as the i-th cell
% 
%     plot(t(1:(end-tbuff(i,1)+1)), mpu_diff_cell{i})
%     xline(130*T);xline(300*T);xline(483*T);xline(605*T); xline(817*T); xline(950*T);
%     xline(1200*T); xline(1331*T);xline(1531*T);xline(1652*T); xline(1838*T); xline(1925*T);
%     title(['Plot of Model subtracted from MPU Trial ', num2str(i)])
%     xlim([0 40]); ylim([-0.5 0.5])
% 
% 
%     subplot(212)  % Matlab diff plots
% % Calculate the difference vector for mat_cln
%     mat_diff = mat_cln(tbuff(i,2):end,i)' - model(1:length(mpu_cln(:,i)) - tbuff(i,2)+1);
%     mat_diff_cell{i} = mat_diff; % Store the vector in the cell array
% 
%     plot(t(1:(end-tbuff(i,2)+1)), mat_diff_cell{i})
%     xline(130*T);xline(300*T);xline(483*T);xline(605*T); xline(817*T); xline(950*T);
%     xline(1200*T); xline(1331*T);xline(1531*T);xline(1652*T); xline(1838*T); xline(1925*T);
%     title(['Plot of Model subtracted from Matlab Trial ', num2str(i)])
%     xlim([0 40]); ylim([-0.5 0.5])
% end

%%
 %Plot of Difference of Model & Data

 figure;
 subplot(211) % MPU Diff plots
 i = 9;
 mpu_diff =  mpu_cln(tbuff(i,1):end,i)' - model(1:length(mpu_cln(:,i)) - tbuff(i,1)+1);
 mpu_diff_cell{i} = mpu_diff; % Store the vector as the i-th cell

 plot(t(1:(end-tbuff(i,1)+1)), mpu_diff_cell{i})

 xline(2.06);xline(6.1);xline(9);xline(12.16); xline(15.3); xline(16.76);
 xline(22.5); xline(22.5);xline(26.7);xline(29.8); xline(33.4); xline(1800*T); xline(1875*T);

 title(['Plot of Model subtracted from MPU Trial ', num2str(i)])
 title(['Plot of Model subtracted from MPU Trial'])
 xlim([0 40]); ylim([-0.5 0.5])
xlabel('Time (s)')
ylabel('Change in acceleration (m/s^2)')
grid on; box


% figure;
 subplot(212)  % Matlab diff plots
 % Calculate the difference vector for mat_cln
 i = 8;
 mat_diff = mat_cln(tbuff(i,2):end,i)' - model(1:length(mpu_cln(:,i)) - tbuff(i,2)+1);
 mat_diff_cell{i} = mat_diff; % Store the vector in the cell array

 plot(t(1:(end-tbuff(i,2)+1)), mat_diff_cell{i})
 xline(1.8);xline(300*T);xline(8.74);xline(605*T); xline(15.44); xline(16.52);
 xline(23.26); xline(1331*T);xline(30.04);xline(1652*T);xline(36.16); xline(1838*T);
 title(['Plot of Model subtracted from Matlab Trial ', num2str(i)])
 title(['Plot of Model subtracted from Matlab Trial'])
 xlim([0 40]); ylim([-0.5 0.5])
xlabel('Time (s)')
ylabel('Change in acceleration (m/s^2)')
grid on; box



%% Determing variance of each stage for each trial
% Define offsets
offsets = [100, 150, 200, 100, 225, 110, 265, 125, 225, 150, 150, 75]; 

% Initialize the result matrices
var_mpu = zeros(13, 15);  % 13 rows for the 13 intervals, 15 columns for the 15 trials
var_mat = zeros(13, 15);

% Loop over each trial
for i = 1:length(tbuff)
    % Initialize t_mpu and t_mat matrices for the current trial
    t_mpu(1, i) = tbuff(i, 1) + offsets(1);  % First time point: tbuff(i, 1) + offsets(1)
    t_mat(1, i) = tbuff(i, 2) + offsets(1);  % First time point: tbuff(i, 2) + offsets(1)
    
    % Calculate subsequent time points by adding cumulative offsets
    for j = 2:12
        t_mpu(j, i) = t_mpu(j-1, i) + offsets(j);  % Add the offset to the previous time point
        t_mat(j, i) = t_mat(j-1, i) + offsets(j);  % Add the offset to the previous time point
    end
    
    % Loop over the 13 time intervals for the current trial
    for j = 1:13
        % If it's the first interval, calculate variance from tbuff(i,1) to t_mpu(1, i)
        if j == 1
            var_mpu(j, i) = var(mpu_cln(tbuff(i, 1):t_mpu(j, i))); 
            var_mat(j, i) = var(mat_cln(tbuff(i, 2):t_mat(j, i))); 
        % If it's the last interval, calculate variance from t_mpu(12, i) to the end of the data
        elseif j == 13
            var_mpu(j, i) = var(mpu_cln(t_mpu(j-1, i):end)); 
            var_mat(j, i) = var(mat_cln(t_mat(j-1, i):end)); 
        % For intermediate intervals, calculate variance from previous t_mpu to current t_mpu
        else
            var_mpu(j, i) = var(mpu_cln(t_mpu(j-1, i) + 1:t_mpu(j, i))); 
            var_mat(j, i) = var(mat_cln(t_mat(j-1, i) + 1:t_mat(j, i))); 
        end
    end
end

%%
tmpu = [1613, 1664, 1620, 1605, 1616, 1635, 1621, 1622, 1581, 1581, 1571, 1602, 1595, 1608, 1570];
tmat = [1603, 1633, 1615, 1620, 1625, 1605, 1600, 1614, 1670, 1583, 1600, 1574, 1605, 1590, 1592];

for n = 1:15
deltaHmpu(n) = mpu_p(tmpu(n),n);
deltaHmat(n) = mat_p(tmat(n),n);
end


% % Indices to exclude
% exclude_indices = [4, 10, 12, 15];
% 
% % Extract all values except the excluded indices
% result = deltaHmpu;
% result(exclude_indices) = [];
% 
% delVarmpu = var(result)
% delStdmpu = std(result)
% 
% 
% deltaHmat(1) = [];

delhVarmpu = var(deltaHmpu,1)
delhVarmat = var(deltaHmat,1)

delHSTDmpu = std(deltaHmpu)
delHSTDmat = std(deltaHmat)

%%
figure; 
hold on; % Hold to plot all trials on the same axes

for i = 1:15 % Loop over trials
    plot(var_mpu(:,i), 'o-', 'DisplayName', sprintf('Trial %d', i)); % Add trial to the plot
end
ylim([0 0.06])
title('Variance at Each Stage for MPU Trials');
xlabel('Stage');
ylabel('Variance');
legend show; % Automatically display legend with trial labels

figure; 
hold on; % Hold to plot all trials on the same axes

for i = 1:15 % Loop over trials
    plot(var_mat(:,i), 'o-', 'DisplayName', sprintf('Trial %d', i)); % Add trial to the plot
end
ylim([0 0.06])
title('Variance at Each Stage for MAT Trials');
xlabel('Stage');
ylabel('Variance');
legend show; % Automatically display legend with trial labels
%%
var_mpu_mean = mean(var_mpu,2);
var_mat_mean = mean(var_mat,2);

rowNames = cell(13, 1); % Preallocate a cell array
for i = 1:13
    rowNames{i} = ['Stage ', num2str(i)];
end

var_mean = table(var_mpu_mean*10^2, var_mat_mean*10^2, 'VariableNames', {'MPU', 'Matlab'}, 'RowNames',rowNames)
%%
t1 = 1; t2 = 2.5*50+1;
t3 = 12.2*50+1; t4=24*50+1;
t5 = 33.8*50+1; t6 = 40*50+1;

stage1 = clean_data(t1:t2);
stage2 = clean_data(t2:t3);
stage3 = clean_data(t3:t4);
stage4 = clean_data(t4:t5);
stage5 = clean_data(t5:t6);


var_s1 = var(stage1);
var_s2 = var(stage2);
var_s3 = var(stage3);
var_s4 = var(stage4);
var_s5 = var(stage5);

%%
figure; 
hold on; % Hold to plot all trials on the same axes

for i = 1:15 % Loop over trials
    plot(t, mpu_p(:,i), 'DisplayName', sprintf('Trial %d', i)); % Add trial to the plot
end
% ylim([0 0.06])
title('Height for different MPU Trials');
xlabel('Time');
ylabel('Height (m)');
legend show; % Automatically display legend with trial labels
grid on;
%%
figure; 
hold on; % Hold to plot all trials on the same axes

for i = 1:15 % Loop over trials
    plot(t, mat_p(:,i), 'DisplayName', sprintf('Trial %d', i)); % Add trial to the plot
end
% ylim([0 0.06])
title('Height for different MAT Trials');
xlabel('Time');
ylabel('Height (m)');
legend show; % Automatically display legend with trial labels
grid on;
toc