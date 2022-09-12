% MEC 2015 System identification
% Prof Helon Ayala

% 04_RLS - recursive least squares
% application to time-varying systems

clc
clear
close all
rng(1) % allows reproduction
% - code begins

%% generate input/output data

% number of measurements
N = 1250;
t = (1:N)';
sig = 1e-2; % std of noise corrupting measurements

% model parameters
a1 = -1.51*ones(N,1);
a2 = 0.55*ones(N,1);
b1 = 0.6*ones(N,1);
b2 = 0.16*ones(N,1);

% add changes in the parameters
a1(251:end)  = 0.8*a1(1);
a2(501:end)  = 1.2*a2(1);
b1(751:end)  = 0.3*b1(1);
b2(1001:end) = 2*b2(1);

% plot parameters
figure
subplot(411),plot(t,a1,'-k','Linewidth',2), title('a1')
subplot(412),plot(t,a2,'-k','Linewidth',2), title('a2')
subplot(413),plot(t,b1,'-k','Linewidth',2), title('b1')
subplot(414),plot(t,b2,'-k','Linewidth',2), title('b2')

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
legend('y','y (noise)')
subplot(2,1,2)
plot(t,u)
legend('u')
xlabel('Time')

%% perform RLS

% initialize parameters
P  = 1e4*eye(4);  % P initial condition
th = 0*ones(4,1); % theta initial condition
lambda = 0.98;     % forgetting factor
yh = zeros(N,1);  % prediction
e  = zeros(N,1);  % prediction error

a1h = zeros(N,1);
a2h = zeros(N,1);
b1h = zeros(N,1);
b2h = zeros(N,1);

% strategy = 1; % RLS
% strategy = 2; % random walk
strategy = 3; % forgetting factor

for k=3:N
    % Step 2
    phi = [-y(k-1) -y(k-2) u(k-1) u(k-2)]'; % regression vector
    
    % Step 3
    yh(k) = phi'*th;
    e(k) = y(k) - yh(k);
    
    % Step 4
    switch strategy
        case 3   % forgetting factor
            K = P*phi/(lambda+phi'*P*phi);
        otherwise
            K = P*phi/(1+phi'*P*phi);
    end
    
    % Step 5
    th = th + K * e(k);
    
    % Step 6    
    switch strategy
        case 1 % RLS
            P = P-K*(P*phi)'; 
        case 2 % random walk
            P = P-K*(P*phi)'; 
            p = length(th)+1;
            q = trace(P)/p;
            Q = q*eye(size(P));
            P = P + Q;
        case 3 % forgetting factor
            P = 1/lambda*(P-K*(P*phi)');
    end
    
    % store estimated parameters for plotting
    a1h(k) = th(1);
    a2h(k) = th(2);
    b1h(k) = th(3);
    b2h(k) = th(4);
end

%% final plots

% plot output: real vs estimated and residual
figure
subplot(2,1,1)
plot(t,[yr yh])
legend('y (real)','y (estimated)')
subplot(2,1,2)
plot(t,e)
legend('error')

% plot output: real vs estimated
figure
subplot(411),plot(t,a1,'-k','Linewidth',2), hold on, plot(t, a1h,'.-k'), title('a1'), legend('real','estimated');
subplot(412),plot(t,a2,'-k','Linewidth',2), hold on, plot(t, a2h,'.-k'), title('a2'), legend('real','estimated');
subplot(413),plot(t,b1,'-k','Linewidth',2), hold on, plot(t, b1h,'.-k'), title('b1'), legend('real','estimated');
subplot(414),plot(t,b2,'-k','Linewidth',2), hold on, plot(t, b2h,'.-k'), title('b2'), legend('real','estimated');
