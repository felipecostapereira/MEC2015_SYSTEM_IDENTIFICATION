% based on sysid.m file from casadi example pack
% updated by Lucas Souza and Helon Ayala 06-2021
% Description:
% in this file we proceed the grey-box identification of state-space
% nonlinear models.
% it is possible to adapt it to different cases, by loading different data
% and changing the state equations as below

clc
clear
close all

addpath(genpath('casadi/'))

import casadi.*

%%%%%%%%%%%% load data %%%%%%%%%%%%%%%%%%%%%
% File to load
load('DATA_EMPS')
% Variables are:
% qm = motor position (measured through the motor encoder)
% qg = the reference position
% t = time
% vir = motor voltage (output of the controller)

% Construction of the vector measurements
Force1 = gtau*vir; % Motor force

decimate = 5; % see Forgione, Piga, EJC 59 (2021) 69-81
u_data = Force1(1:decimate:end);
y_data = qm(1:decimate:end);
t      = t(1:decimate:end); 

% plot(t,y_data)
%%
N  = length(t);  % Number of samples
Ts = round(t(2)-t(1),3);  % sampling time (seconds)
fs = 1/Ts;       % Sampling frequency [hz]

x0 = DM([qm(1),0]); % initial condition for simulation

%%%%%%%%%%%% MODELING %%%%%%%%%%%%%%%%%%%%%
q  = MX.sym('q');
dq = MX.sym('dq');
u  = MX.sym('u');

states = [q;dq];
controls = u;

Mn    = MX.sym('Mn');
Fvn   = MX.sym('Fvn');
Fcn   = MX.sym('Fcn');
ofstn = MX.sym('ofstn');

params   = [Mn;Fvn;Fcn;ofstn];
parammax = [150; 300; 40;  15];
parammin = [ 30; 100;  0; -15];

nparam = length(params);
param_guess = rand(nparam,1);

lbx = zeros(nparam,1);
ubx = ones(nparam,1);

M    = denorm(Mn,   parammax(1),parammin(1));
Fv   = denorm(Fvn,  parammax(2),parammin(2));
Fc   = denorm(Fcn,  parammax(3),parammin(3));
ofst = denorm(ofstn,parammax(4),parammin(4));

rhs = [dq; (u-Fv*dq-Fc*sign(dq)-ofst)/M];

% Form an ode function
ode = Function('ode',{states,controls,params},{rhs});

%%%%%%%%%%%% Creating a simulator %%%%%%%%%%
N_steps_per_sample = 10;
dt = 1/fs/N_steps_per_sample;

% Build an integrator for this system: Runge Kutta 4 integrator
k1 = ode(states,controls,params);
k2 = ode(states+dt/2.0*k1,controls,params);
k3 = ode(states+dt/2.0*k2,controls,params);
k4 = ode(states+dt*k3,controls,params);

states_final = states+dt/6.0*(k1+2*k2+2*k3+k4);

% Create a function that simulates one step propagation in a sample
one_step = Function('one_step',{states, controls, params},{states_final});

X = states;
for i=1:N_steps_per_sample
    X = one_step(X, controls, params);
end
%
% % Create a function that simulates all step propagation on a sample
one_sample = Function('one_sample',{states, controls, params}, {X});
%
% speedup trick: expand into scalar operations
one_sample = one_sample.expand();

%%%%%%%%%%%% Simulating the system %%%%%%%%%%

all_samples = one_sample.mapaccum('all_samples', N);

%%%%%%%%%%%% Identifying the simulated system %%%%%%%%%%
opts = struct;
% opts.ipopt.max_iter = 15;
% opts.ipopt.print_level = 3;%0,3
% opts.print_time = 1;
opts.ipopt.acceptable_tol = 1e-4;
opts.ipopt.acceptable_obj_change_tol = 1e-4;

single_multiple = 0;
switch single_multiple
    case 1           
        %%%%%%%%%%%% single shooting strategy %%%%%%%%%%
        X_symbolic = all_samples(x0, u_data, repmat(params,1,N));
        
        e = y_data-X_symbolic(1,:)';
        
        J = 1/N*dot(e,e);
        nlp = struct('x', params, 'f', J);
        
        solver = nlpsol('solver', 'ipopt', nlp, opts);
        
        sol = solver('x0', param_guess, 'lbx', lbx,'ubx', ubx);
        
        % parametros identificados:
        paramhat = sol.x.full;
    otherwise
        %%%%%%%%%%%% multiple shooting strategy %%%%%%%%%%
        % % All states become decision variables
        X = MX.sym('X', 2, N);

        res = one_sample.map(N, 'thread', 4);
        Xn = res(X, u_data', repmat(params,1,N));

        gaps = Xn(:,1:end-1)-X(:,2:end);

        e = y_data-Xn(1,:)';

        V = veccat(params, X);

        J = 1/N*dot(e,e);

        nlp = struct('x',V, 'f',J,'g',vec(gaps));

        % Multipleshooting allows for careful initialization
        yd = diff(y_data)*fs;
        X_guess = [ y_data  [yd;yd(end)]]';

        param_guess = [param_guess(:);X_guess(:)];

        solver = nlpsol('solver','ipopt', nlp, opts);

        sol = solver('x0',param_guess,'lbg',0,'ubg',0);
        solx = sol.x.full;
        paramhat = solx(1:nparam);
end

%% analisa resultado

Mhat    = denorm(paramhat(1),parammax(1),parammin(1));
Fvhat   = denorm(paramhat(2),parammax(2),parammin(2));
Fchat   = denorm(paramhat(3),parammax(3),parammin(3));
ofsthat = denorm(paramhat(4),parammax(4),parammin(4));

disp('Parametros identificados:')
[Mhat, Fvhat, Fchat, ofsthat]

disp('Parametros IDIM:')
paramhatIDIM = [95.1089, 203.5034, 20.3935, -3.1648]';
paramhatIDIM' 
paramhatIDIM = normalize(paramhatIDIM,parammax,parammin); % simulation


%% compare both solutions (IDIM vs. CASADI)

Xhat = all_samples(x0, u_data, repmat(paramhat,1,N));
Xhat = Xhat.full;
yhat1 = Xhat(1,:)';

XhatIDIM = all_samples(x0, u_data, repmat(paramhatIDIM,1,N));
XhatIDIM = XhatIDIM.full;
yhatIDIM1 = XhatIDIM(1,:)';

% saving vars to plot
save('r_coulomb.mat','t','y_data','yhat1','yhatIDIM1')

figure
hold on
plot(t,y_data,'k-','linewidth',1.5)
plot(t,yhat1,'r--','linewidth',1.5)
plot(t,yhatIDIM1,'g--','linewidth',1.5)
grid on
xlabel('time')
legend({'real','casadi','IDIM'},'location','best')

figure
hold on
plot(t,y_data-yhat1,'r-','linewidth',1.5)
plot(t,y_data-yhatIDIM1,'g--','linewidth',1.5)
grid on
xlabel('time')
ylabel('error')
legend({'casadi','IDIM'},'location','best')
