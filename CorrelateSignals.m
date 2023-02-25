%% Code written by Amit Chowdhury
%% Indian Institute of Technology, Jammu
%% email: 2022psp0004@iitjammu.ac.in
%% Date: 25th Feb 2023

clc
clear
close all

% Define the parameters of the signals
Fs = 1000;       % Sampling frequency
T = 1/Fs;        % Sampling period
n = 1000;        % Length of signals
t = (0:n-1)*T;   % Time vector

% Generate random signals
x = randn(n, 1); % signal 1
y = randn(n, 1); % signal 2

% Compute correlation coefficient
rho = 0.8; % desired correlation coefficient
C = [1, rho; rho, 1]; % correlation matrix
R = chol(C); % Cholesky factorization

% Generate correlated signals
z = R * [x, y].'; % correlated signals
x_corr = z(1, :).'; % correlated signal 1
y_corr = z(2, :).'; % correlated signal 2

figure;
subplot(2,2,1);            % plots of uncorrelated signals
plot(t, z(1,:), 'b');
xlabel('Time (s)');
ylabel('Signal 1');
subplot(2,2,2);
plot(t, z(2,:), 'r');
xlabel('Time (s)');
ylabel('Signal 2');

subplot(2,2,3);          % plots of correlated signals
plot(t, x_corr, 'b');
xlabel('Time (s)');
ylabel('Signal 1');
subplot(2,2,4);
plot(t, y_corr, 'r');
xlabel('Time (s)');
ylabel('Signal 2');

figure;
subplot(2,1,1);
plot(xcorr(x(:),y(:)));   %% xcorr(x,y) returns the cross-correlation of two discrete-time sequences.
title('CrossCorr betw two uncorrelated signals');
subplot(2,1,2);
plot(xcorr(z(1,:),z(2,:)));
title('CrossCorr betw two correlated signals');


% In this example, we first generate two random signals x and y using the randn function. We then compute the correlation 
%coefficient rho that we want to achieve between the two signals. We create the correlation matrix C using rho, and then 
%compute its Cholesky factorization R using the chol function.

% Finally, we generate the correlated signals x_corr and y_corr by multiplying the Cholesky factor R by a matrix containing 
%the original signals x and y. The resulting matrix z contains the correlated signals, and we extract each signal by selecting 
%the appropriate row.

% Note that you can adjust the value of rho to control the strength of the correlation between the two signals. 
%A value of rho = 1 corresponds to perfect correlation, while a value of rho = 0 corresponds to no correlation.
