clc
clear
close all


load('r_coulomb.mat')
load('r_tustin.mat')

figure
hold on
plot(t(1:50:end),y_data(1:50:end),'go','linewidth',2, 'MarkerSize', 8)
plot(t,yhat1,'r--','linewidth',1.5)
plot(t,yhatIDIM1,'b--','linewidth',1.5)
plot(t,yhat2,'r-','linewidth',1.5)
plot(t,yhatIDIM2,'b-','linewidth',1.5)
grid on
xlabel('time')
ylabel('y - position')
legend({'real','Coulomb(casadi)','Coulomb(IDIM)','Tustin(casadi)','Tustin(IDIM)',},'location','best')

figure
hold on
plot(t,y_data-yhat1,'r--','linewidth',1.5)
plot(t,y_data-yhatIDIM1,'b--','linewidth',1.5)
plot(t,y_data-yhat2,'r-','linewidth',1.5)
plot(t,y_data-yhatIDIM2,'b-','linewidth',1.5)
grid on
xlabel('time')
ylabel('error')
legend({'Coulomb(casadi)','Coulomb(IDIM)','Tustin(casadi)','Tustin(IDIM)',},'location','best')
