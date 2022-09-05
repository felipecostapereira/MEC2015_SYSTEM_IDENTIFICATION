% MEC 2015 System identification
% Prof Helon Ayala

% 01_MQB

clc
clear
close all
rng(1) % permite reproducibilidade

load utra % training data created (multisine)

% ---- code begins

%% - generate artificial data
% spring-mass-damper system
th = [-1.895 0.9048 0.0004833 0.0004675]'; % true vector of parameters

% simulation parameters
N = 1000; % simulation steps
t = (1:N)'; % time vector
sig = 0.02; % noise corrupting measurements

% ESTIMATION
u = utra;
y = zeros(N,1);
for k=3:N
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
uamp = 1;
yr = 0.05;
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

% free run simulation
yhat_TRA_FR = zeros(N,1);
yhat_VAL_FR = zeros(N,1);
yhat_TRA_FR(1:2) = YTRA(1:2); % initialize with measurements
yhat_VAL_FR(1:2) = YVAL(1:2); % initialize with measurements
for k=3:N
    yhat_TRA_FR(k) = -th(1)*yhat_TRA_FR(k-1) + -th(2)*yhat_TRA_FR(k-2) + th(3)*UTRA(k-1) + th(4)*UTRA(k-2);
    yhat_VAL_FR(k) = -th(1)*yhat_VAL_FR(k-1) + -th(2)*yhat_VAL_FR(k-2) + th(3)*UVAL(k-1) + th(4)*UVAL(k-2);
end
yhat_TRA_FR = yhat_TRA_FR(3:end); % remove measurements
yhat_VAL_FR = yhat_VAL_FR(3:end); % remove measurements
%% plot predictions
figure
plot([Y1 yhat_TRA_OSA yhat_TRA_FR])
title('Estimation phase')
legend('Measured','One-step-ahead prediction','Free-run simulation')

figure
plot([Y2 yhat_VAL_OSA yhat_VAL_FR])
title('Validation phase')
legend('Measured','One-step-ahead prediction','Free-run simulation')

%% plot residuals
figure
plot([Y1-yhat_TRA_OSA Y1-yhat_TRA_FR])
title('Residual in estimation phase')
legend('One-step-ahead prediction','Free-run simulation')

figure
plot([Y2-yhat_VAL_OSA Y2-yhat_VAL_FR])
title('Residual in validation phase')
legend('One-step-ahead prediction','Free-run simulation')


%% calculate metrics
R2_TRA_OSA = mult_corr(Y1,yhat_TRA_OSA);
R2_VAL_OSA = mult_corr(Y2,yhat_VAL_OSA);
R2_TRA_FR  = mult_corr(Y1,yhat_TRA_FR);
R2_VAL_FR  = mult_corr(Y2,yhat_VAL_FR);

fprintf('R2 obtained (1 is perfect) \n')
fprintf('R2_TRA_OSA \t R2_VAL_OSA \t R2_TRA_FR \t R2_VAL_FR\n')
fprintf('%0.4f \t\t %0.4f \t\t %0.4f \t %0.4f\n',R2_TRA_OSA,R2_VAL_OSA,R2_TRA_FR,R2_VAL_FR)


%% correlation based tests
plot_xcorrel(Y1-yhat_TRA_OSA,UTRA(3:end))
