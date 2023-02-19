%%% MATLAB Code for Direction of Arrival in aerial. 
%%% Code written by MohammadReza Jabbari-Email: Mo.re.jabbari@gmail.com
%%% Department of Electrical and Computer Engineering
%%% Isfahan University of Technology-IRAN

clc
clear
close all

format compact    

DOA = [ 60  -10];      %Direction of arrival (Degree)
T   = 200;         %Snapshots (or Samples)
K   = length(DOA); %The number of signal source(or traget)
Nr  = 6;          %Number of sensors
lambda = 1;      %Wavelength
d   = lambda/2;    %Receiver's antennas spacing
SNR = 10;           %Signal to Noise Ratio (dB)
A = zeros(Nr,K);   %Steering Matrix

for k=1:K
    A(:,k) = exp(-1j*2*pi*d*sind(DOA(k))*(0:Nr-1)'/lambda); %Assignment matrix
end

Vj = diag(sqrt((10.^(SNR/10))/2));
s = Vj* ( randn(K,T) + 1j*randn(K,T) );
noise = sqrt(1/2)*(randn(Nr,T)+1j*randn(Nr,T));

X = A*s;
X = X+noise;      %Insert Additive White Gaussain Noise (AWGN)

%% MUSIC (MUltiple SIgnal Classification)

Rx = cov(X');                     %Data covarivance matrix

[eigenVec,eigenVal] = eig(Rx);    %Find the eigenvalues and eigenvectors of Rx

Vn = eigenVec(:,1:Nr-K);          %Estimate noise subspace (Note that eigenvalues sorted ascendig on columns of "eigenVal")

theta = -90:0.5:90;       %Grid points of Peak Search

for i=1:length(theta)

    SS = zeros(Nr,1);

    SS = exp(-1j*2*pi*d*(0:Nr-1)'*sind(theta(i))/lambda);

    PP = SS'*(Vn*Vn')*SS;

    Pmusic(i) = 1/ PP;

end

Pmusic=abs(Pmusic./max(Pmusic));
Pmusic1 = real(10*log10(Pmusic)); %Spatial Spectrum function

[pks, locs] = findpeaks(Pmusic,theta,'SortStr','descend','Annotate','extents');

MUSIC_Estim = sort(locs(1:K))

figure;
plot(theta,Pmusic1,'-b',locs(1:K), pks(1:K),'r*'); hold on
text(locs(1:K)+2*sign(locs(1:K)), pks(1:K),num2str(locs(1:K)'))
xlabel('Angle \theta (degree)'); ylabel('Spatial Power Spectrum P(\theta) (dB)')
title('DOA estimation based on MUSIC algorithm ')
xlim([min(theta) max(theta)])
grid on
