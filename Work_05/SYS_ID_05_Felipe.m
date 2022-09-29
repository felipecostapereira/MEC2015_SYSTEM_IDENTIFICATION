clc 
clear 
close all

addpath(genpath('Contsid'))

% time R2 R3 An2 An3 DC_Torque
rawData = readmatrix('measured_data.csv');

Ts = 0.01;
fs = 1/Ts;
[resData, t] = resample(rawData(:,2:end),rawData(:,1),fs);

figure
subplot(2,1,1)
plot(t,resData(:,5))
ylabel('Input - DC Torque')
subplot(2,1,2)
plot(t,resData(:,2))
ylabel('Output - velocity of J2')

% iddata 
z = iddata(resData(:,1),resData(:,5), Ts, 'interSample', 'zoh');
[B,A] = butter(2,.1);
zf = idfilt(z,{B,A});

u = z.InputData;
uf = zf.InputData;
yf = zf.OutputData;
y = z.OutputData;

% plotting filtered data
figure
subplot(2,1,1)
plot(t,[u, uf])
ylabel('Input')
legend ('u', 'u filtered')
subplot(2,1,2)
plot(t,[u, uf])
ylabel('Output')
legend ('y', 'y filtered')

% plotting u and y power spectrum
figure
subplot(2,1,1)
pspectrum(u, fs)
hold on
pspectrum(uf, fs)
ylabel('Input')
legend ('u', 'u filtered')
subplot(2,1,2)
pspectrum(y, fs)
hold on
pspectrum(yf, fs)
hold on

% Looping over params to proceed system identification
% Measuring quality of fit using R2 Score
best_r2 = -100;
for np = 2:8
    for nz = 0:np-1
        Ghat = tfrivc(zf,np,nz);
        ye = lsim(Ghat, uf, t);  
        r2 = mult_corr(yf,ye);            
        fprintf('np=%d, nz=%d, R2= %0.4f ', np, nz, r2);
        if r2 > best_r2                
            fprintf('(best so far)\n');
            best_r2 = r2;
            best_Ghat = Ghat;                
        else
            fprintf('\n');
        end        
    end
end
fprintf('>>>>>>>>>>>>>>>>>>> Best R2=%.4f\n',best_r2);

% simulating system (best Ghat)
ye = lsim(best_Ghat, uf, t);       
ze = iddata(ye, uf, Ts);

% measured vs estimated
figure
plot(t,[yf, ye]);
ylabel('y');
legend('y (filtered)', 'y (estimted)')

figure
pspectrum(y, fs)
hold on
pspectrum(yf, fs)
hold on
pspectrum(ye, fs)
ylabel('Output')
legend ('y', 'y filtered', 'y estimated')


