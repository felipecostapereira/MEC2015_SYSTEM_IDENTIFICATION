clc 
clear 
close all

% time R2 R3 An2 An3 DC_Torque
rawData = dlmread('measured_data.csv');

% plot(rawData(:,1),rawData(:,5))
% plot(diff(rawData(:,1)))

Ts = 1e-2;
fs = 1/Ts;
[resData, t] = resample(rawData(:,2:end),rawData(:,1),fs);

plot(diff(t))

figure
subplot(2,1,1)
plot(t,resData(:,5))
subplot(2,1,2)
plot(t,resData(:,2))

% ind = t > 10 & t < 40;

z = iddata(resData(:,1),resData(:,5));
[B,A] = butter(2,1e-1);
zf = idfilt(z,{B,A});

figure
plot(z,zf)

