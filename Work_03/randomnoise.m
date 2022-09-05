
close all
clc


% filtered random noise 
Ndata = 2048;
fsample = 1000;
Cutoff = 0.1; % cutoff normalized freq
[b,a] = butter(6,Cutoff*2);
u = filter(b,a,randn(Ndata,1));

Um=fft(u)/sqrt(Ndata);       % spectrum of the actual generate multisine
LinesPlot=1:floor(Ndata/2);  % lines to be plotted
f=(0:Ndata-1)*fsample/Ndata;  % frequencies

figure
subplot(1,2,2)
plot(f(LinesPlot),db(Um(LinesPlot)),'k')
xlabel('Frequency (Hz)'),ylabel('Amplitude (dB)')
subplot(1,2,1)
plot((0:(Ndata-1))'/fsample,u,'k')
xlabel('time [s]'),ylabel('u(t)')
