% MEC 2015 System identification
% Prof Helon Ayala

% 04_RLS - recursive least squares

clc
clear
close all
rng(1) % allows reproduction
% - code begins

%% generate input/output data

% number of measurements
N = 200;
t = (1:N)';
sig = 0*1e-2; % std of noise corrupting measurements

% model parameters
a1 = -1.51*ones(N,1); 
a2 = 0.55*ones(N,1);
b1 = 0.6*ones(N,1);
b2 = 0.16*ones(N,1);

% init
yr = zeros(N,1); % model (no noise)
u = zeros(N,1);
for i=1:N
    if rand>0.5
        u(i) = 1;
    else
        u(i) = -1;
    end
end

for k=3:N
    yr(k) = - a1(k)*yr(k-1) - a2(k)*yr(k-2) + b1(k)*u(k-1) + b2(k)*u(k-2);
end

y = yr + sig.*randn(N,1); % measurements (corrupted by noise)

% plot input/output data
figure
subplot(2,1,1)
plot(t,[yr y])
legend('y')
subplot(2,1,2)
plot(t,u)
legend('u')
xlabel('Time')

%% perform RLS

% initialize parameters
P = 1e0*eye(4);
th = -100*ones(4,1);
yh = zeros(N,1); % prediction
e = zeros(N,1);  % prediction error
trP = zeros(N,1);  % prediction error
trP(1) = trace(P);
trP(2) = trace(P);

a1h = zeros(N,1);
a2h = zeros(N,1);
b1h = zeros(N,1);
b2h = zeros(N,1);

for k=3:N
    % Step 2
    phi = [-y(k-1) -y(k-2) u(k-1) u(k-2)]'; % regression vector
    
    % Step 3
    yh(k) = phi'*th;
    e(k) = y(k) - yh(k);
    
    % Step 4
    K = P*phi/(1+phi'*P*phi);
    
    % Step 5
    th = th + K * e(k);
    
    % Step 6
    P = P-K*(P*phi)';
    
    % store estimated parameters for plotting
    a1h(k) = th(1);
    a2h(k) = th(2);
    b1h(k) = th(3);
    b2h(k) = th(4);
    trP(k) = trace(P);
end

%% final plots

% plot output: real vs estimated and residual
figure
plot(t,[yr yh e])
legend('y (real)','y (estimated)','error')

% plot covariance matrix trace
figure
plot(t,trP)

% plot output: real vs estimated
figure
subplot(411),plot(t,a1,'-k','Linewidth',2), hold on, plot(t, a1h,'.-k'), title('a1'), legend('real','estimated');
subplot(412),plot(t,a2,'-k','Linewidth',2), hold on, plot(t, a2h,'.-k'), title('a2'), legend('real','estimated');
subplot(413),plot(t,b1,'-k','Linewidth',2), hold on, plot(t, b1h,'.-k'), title('b1'), legend('real','estimated');
subplot(414),plot(t,b2,'-k','Linewidth',2), hold on, plot(t, b2h,'.-k'), title('b2'), legend('real','estimated');
