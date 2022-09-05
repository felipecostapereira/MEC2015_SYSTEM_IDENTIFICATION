% MEC 2015 System identification
% Prof Helon Ayala

% 01_MQB

clc
clear
close all
rng(1) % permite reproducibilidade

% ---- code begins

%% - generate artificial data
% spring-mass-damper system

m = 0.1; b= 1; k = 10; % parameters
s = tf('s');
G = 1/(m*s^2+b*s+k); % continuous model
Ts = 0.01; % sampling time
Gd = c2d(G,Ts) % discrete model

th = [-1.895 0.9048 0.0004833 0.0004675]'; % true vector of parameters

% simulation parameters
N = 1000; % simulation steps
t = (1:N)'; % time vector
sig = 0.01; % noise corrupting measurements

% TRAINING
uamp = 10;
yr = 1;
u = uamp*ones(N,1);
y = zeros(N,1);
for k=3:N
    if y(k-1) >=  yr
        u(k) = -uamp;
    elseif y(k-1) <=  -yr
        u(k) = uamp;
    else
        u(k) = u(k-1);        
    end
    y(k) = -th(1)*y(k-1) + -th(2)*y(k-2) + th(3)*u(k-1) + th(4)*u(k-2);
end
UTRA = u;
YTRA = y + sig*randn(N,1);

% plot training data
figure
subplot(2,1,1)
plot(t,UTRA)
title('input')
legend('u');
subplot(2,1,2)
plot(t,[YTRA y])
title('output')
legend('y','y measured (+noise)');

% VALIDATION
u = uamp/4*sin(2*pi*t*Ts) + uamp/4*sin(pi/2*t*Ts) + uamp/4*sin(pi*t*Ts) + uamp/4*sin(pi/4*t*Ts);
y = zeros(N,1);
for k=3:N
    y(k) = -th(1)*y(k-1) + -th(2)*y(k-2) + th(3)*u(k-1) + th(4)*u(k-2);
end
UVAL = u;
YVAL = y + sig*randn(N,1);

% plot validation data
figure
subplot(2,1,1)
plot(t,UVAL)
title('input')
legend('u');
subplot(2,1,2)
plot(t,[YVAL y])
title('output')
legend('y','y measured (+noise)');

%% - perform identification

% construct regression matrix
%        [  -y(k-1)            -y(k-2)        u(k-1)        u(k-2)     ]
Phi    = [-YTRA(2:end-1)  -YTRA(1:end-2)  UTRA(2:end-1)   UTRA(1:end-2)];
PhiVAL = [-YVAL(2:end-1)  -YVAL(1:end-2)  UVAL(2:end-1)   UVAL(1:end-2)]; % only for validation
Y1 = YTRA(3:end); % targets
Y2 = YVAL(3:end); % targets

% estimate parameters:
th_hat = (Phi'*Phi)^(-1)*Phi'*Y1; % batch least squares

% one step ahead prediction
yhat_TRA_OSA = Phi*th_hat;
yhat_VAL_OSA = PhiVAL*th_hat;

%% plot predictions
figure
plot([Y1 yhat_TRA_OSA])
title('Estimation phase')
legend('Measured','One-step-ahead prediction')

figure
plot([Y2 yhat_VAL_OSA])
title('Validation phase')
legend('Measured','One-step-ahead prediction')









