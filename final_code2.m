% %% IN THE NAME OF THE MERCIFUL
clear all;
close all;
clc;

Titer=2000;
label = zeros(Titer,1201);
cov_data = cell(Titer,1);
for iter=1:Titer

    %% Underwater Channel
    Sensors=6; propSpeed=1500;channelDepth=200; OperatingFrequency = 10e3;
    % Set the range and resolution for the random number
    path_min = 3;          % minimum value of x
    path_max = 17;         % maximum value of x
    x_res = 1;          % resolution of x
    % Generate the random number
    mu = (path_min + path_max) / 2;
    sigma = (path_max - path_min) / 4;
    %%  ********* this will come in 'for loop'  **********
    numpath = round(normrnd(mu, sigma) / x_res) * x_res;
    isopaths = phased.IsoSpeedUnderwaterPaths('ChannelDepth',channelDepth,'BottomLoss',0.2,...
    'NumPathsSource','Property','NumPaths',numpath,'PropagationSpeed',propSpeed);
    channel = phased.MultipathChannel('OperatingFrequency',OperatingFrequency);
    %% Acoustic Beacon Waveform 1
    prf = 1;
    pulseWidth = 10e-3;
    pulseBandwidth = 1/pulseWidth;
    fs = 2*pulseBandwidth;
    wav = phased.RectangularWaveform('PRF',prf,'PulseWidth',pulseWidth,...
    'SampleRate',fs);
    channel.SampleRate = fs;
% Acoustic Beacon Waveform 2
% Define parameters
% f = 6;                  % Frequency of the cosine wave (in Hz)
% A = 1;                  % Amplitude of the cosine wave
% phi = 0;                % Phase offset of the cosine wave (in radians)
% t_start = 0;           % Start time of the cosine wave (in seconds)
% t_end = 2;             % End time of the cosine wave (in seconds)
% t_step = 0.01;         % Time step (in seconds)
% 
% % Generate time vector
% t = t_start:t_step:t_end;
% 
% % Generate cosine wave
% prf2 = A * cos(2*pi*f*t + phi);
% Plot the cosine wave
%plot(t, x);
%xlabel('Time (s)');
% ylabel('Amplitude');
% title('Cosine Wave');
% grid on;
%prf2 = cos(2*pi*6*t);
%pulseWidth2 = 100e-3;
%pulseBandwidth2 = 1/pulseWidth2;
%fs2 = 2*6;
%wav2 = phased.RectangularWaveform('PRF',prf2,...
%'SampleRate',fs2);
%wav2=prf2;
%channel.SampleRate = fs;
%% Acoustic Beacon 1
% Set the range and resolution for x
    range_x_min = 1000;      % minimum value of x
    range_x_max = 4000;      % maximum value of x
    range_x_res = 5;         % resolution of x
    range_x = range_x_min + (range_x_max - range_x_min) * rand(1);
    range_x = round(range_x / range_x_res) * range_x_res;
% Set the range and resolution for y
    range_y_min = -sqrt(3)*range_x;     % minimum value of y
    range_y_max = sqrt(3)*range_x;      % maximum value of y
    range_y_res = 5;                 % resolution of y
    range_y = range_y_min + (range_y_max - range_y_min) * rand(1);
    range_y = round(range_y/ range_y_res) * range_y_res;
% Set the range and resolution for z
    range_z_min = -1;        % minimum value of z
    range_z_max = -199;      % maximum value of z
    range_z_res = 0.01;      % resolution of z
    range_z = range_z_min + (range_z_max -range_z_min) * rand(1);
    range_z = round(range_z / range_z_res) * range_z_res;
    [azi,el,r] = cart2sph(range_x,range_y,range_z);
% convert radian to degrees
    azi_deg = azi * (180/pi);
% round result to resolution of 0.1
    azi_deg = round(azi_deg * 10)/10 ;
    projector = phased.IsotropicProjector('VoltageResponse',120);
    projRadiator = phased.Radiator('Sensor',projector,...
    'PropagationSpeed',propSpeed,'OperatingFrequency',OperatingFrequency);
    beaconPlat1 = phased.Platform('InitialPosition',[range_x;range_y;range_z],'Velocity',[0; 0; 0]);
%% Passive Towed Array
    hydrophone = phased.IsotropicHydrophone('VoltageSensitivity',-150);
    array = phased.ULA('Element',hydrophone,...
    'NumElements',Sensors,'ElementSpacing',propSpeed/OperatingFrequency/2,...
    'ArrayAxis','y');
    arrayCollector = phased.Collector('Sensor',array,...
    'PropagationSpeed',propSpeed,'OperatingFrequency',OperatingFrequency);
    arrayPlat = phased.Platform('InitialPosition',[0; 0; -20],...
    'Velocity',[0; 0; 0]);%[0; 1; 0]
%Define the receiver amplifier for each hydrophone element.
%Choose a gain of 20 dB and noise figure of 10 dB.
    rx = phased.ReceiverPreamp(...
    'Gain',30,...
    'NoiseFigure',20,...
    'SampleRate',fs,...
    'SeedSource','Property',...
    'Seed',2007);
%% Simulate the Passive Sonar System
    x = wav();
% Set the parameters
    SNR =0;              % SNR in dB
    N = fs;               % number of samples
    sigma = 1;              % standard deviation of noise
% Generate the signal
% signal = ones(N,1);     % actual signal amplitude is 1
% Generate the noise
% noise = sigma*(randn(N,1) + 1i*randn(N,1))/sqrt(2) + ...
    %1i*sigma*(randn(N,1) + 1i*randn(N,1))/sqrt(2) ;  % complex Gaussian noise

% Generate the two independent Gaussian noise sequences
    real_noise = sqrt(sigma/2) * randn(N,1);
    imaginary_noise = sqrt(sigma/2) * randn(N,1);

% Combine the two noise sequences into a complex Gaussian noise sequence
    noise = real_noise + 1i*imaginary_noise;

% Calculate the power of the signal and noise
    signal_power = norm(x)^2/N;
    noise_power = norm(noise)^2/N;
% Calculate the scaling factor for the noise
    scale_factor = sqrt(signal_power/noise_power/10^(SNR/10));
% Scale the noise to achieve the desired SNR
    noise_scaled = scale_factor*noise;
% Add the noise to the signal
    x = x + noise_scaled;
    numTransmits = 6;
    angPassive = zeros(numTransmits,1);
    angAct1 = zeros(numTransmits,1);
    rxsig1 = zeros(size(x,1),Sensors,numTransmits);
    musicspatialspect = phased.MUSICEstimator('SensorArray',array,...
    'PropagationSpeed',propSpeed,'OperatingFrequency',...
    OperatingFrequency,'ScanAngles',-90:0.1:90-1,'DOAOutputPort',true,...
    'NumSignalsSource','Property','NumSignals',1);
% t = (0:length(x)-1)'/fs;
%%
    for i = 1:numTransmits
  % Update array and acoustic beacon positions
      [pos_tx1,vel_tx1] = beaconPlat1(1/prf);
      [pos_rx1,vel_rx1] = arrayPlat(1/prf);
  % Compute paths between the acoustic beacon and array
      [paths1,dop1,aloss1,rcvang1,srcang1] = ...
    isopaths(pos_tx1,pos_rx1,vel_tx1,vel_rx1,1/prf);
      angAct1(i) = rcvang1(1,1);
  % Propagate the acoustic beacon waveform
      tsig1 = projRadiator(x,srcang1);
      rsig1 = channel(tsig1,paths1,dop1,aloss1);
  % Collect the propagated signal
      rsig1 = arrayCollector(rsig1,rcvang1);
 
      RX1=rsig1();
      rxsig1 = rx(RX1);
     [~,angPassive1(i,:)] = musicspatialspect(rxsig1);
    end

    covar=(1/fs).*(rxsig1'*rxsig1);
    covariance(:,:,1)=real(covar);
    covariance(:,:,2)=imag(covar);
    covariance(:,:,3)=angle(covar);

    cov_data{iter}=covariance;
%%
    DOA_ang = round(angPassive1(i,1), 1);
    fprintf("The actual angle of the beacon : %f degrees\n& DOA Estimate using Music : %f degrees\n",azi_deg,DOA_ang);

% Define the range and resolution of the array
    grid = -60:0.1:60;
    if DOA_ang> -1
        index = 601 + DOA_ang/0.1;
    else
        index = 601 - DOA_ang/0.1;
    end
% Set the index to 1
    label(iter,index) = 1;
% Display the array
%disp(label);
end