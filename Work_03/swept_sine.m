% swept sine signal
  
% param
Ts = 1/1e3;
T0 = 1; % period [s]
f0 = 1/T0; % frequency
k1 = 1; % [k1*f0,k2*f0] lowest and highest frequency respectively
k2 = 42; % k2 > k1, natural numbers
a = pi*(k2-k1)*f0^2;
b = 2*pi*k1*f0;
A = 10; % amplitude

tvec= 0:Ts:T0;
% signal
u = A*sin((a*tvec+b).*tvec);
plot(tvec,u)
xlabel('time [s]')
ylabel('u(t)')
