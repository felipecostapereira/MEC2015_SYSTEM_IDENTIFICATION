% MEC 2015 System identification
% Prof Helon Ayala

% 05_CONTSID

clc
clear
close all
rng(1) % permite reproducibilidade

addpath(genpath('Contsid'))

% define model parameters
m = 0.1; b= 0.1; k = 10; % parameters
s = tf('s');
G = 1/(m*s^2+b*s+k); % continuous model

% define excitation signal
Ts = 0.001; % sampling time
fs = 1/Ts;
Tmax = 10; % maximum time
fmax = 10;  % multisine up to 10 Hz
A = 10; % amplitude
ninp = 1; % number of inputs = 1
 
[u,t] = multiSine(fs, fmax, Tmax, A, ninp);

% simulate system
N = length(u);
y = lsim(G,u,t) + 0.1*randn(N,1);

ze = iddata(y,u,Ts,'InterSample','zoh');

figure
plot(ze)

% estimate model using CONTSID
np = 2;
nz = 0;
Ghat = tfrivc(ze,np,nz,'TdMax',0);

% frequency response comparison
figure
bode(G,Ghat)
legend('G','Ghat')

% poles comparison
roots(G.Denominator{1})
roots(Ghat.Denominator)

% model output comparison
ye = lsim(Ghat,u,t);

figure
plot(t,[y,ye]);
legend('measured','simulated')



