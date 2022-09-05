% Chapter 2 Exercise 20
% Generation of a multisine with flat amplitude spectrum
%
% Copyright: 
% Johan Schoukens, Rik Pintelon, and Yves Rolain 
% Vrije Universiteit Brussels, Pleinlaan 2, 1050 Brussels, Belgium
%
% 1 December 2010

clear all
% define the input parameters
fsample=100;Ts=1/fsample;      % sample frequency and sample period
Ndata1=1000;                    % length of the signal
t=[0:Ndata1-1]*Ts;              % time axis
Nsines=50;                     % number of sines
f=[0:Ndata1-1]*fsample/Ndata1;  % multisine frequencies
LinesPlot=[1:floor(Ndata1/2)];  % lines to be plotted

% multisine 2: with a random phase
U2=zeros(Ndata1,1);             % choose random phases
U2(2:Nsines+1)=exp(j*2*pi*rand(Nsines,1));
u2=2*real(ifft(U2));u2=u2/std(u2);
U2m=fft(u2)/sqrt(Ndata1);       % spectrum of the actual generate multisine

plot(t,u2,'k')

utra = u2;

save utra utra