%%% a3Q5
%%% ans to Q no a

clc;
clear;
f=@(t,x) [1.2*x(1)-0.6*x(1)*x(2); -0.8*x(2)+0.3*x(1)*x(2)];
[t, xsol]=ode45(f, [-10 10], [2 1]);


figure;
plot(t,xsol(:,1),'r-O',t,xsol(:,2),'b-*');
title('System of ODE');
xlabel('t');
ylabel('x and y value');
legend('x val','y val');

%%% ans to Q no b

clear;
mu=3; 
f=@(t,y) [y(2);mu*(1-y(1)^2)*y(2)-y(1)];
[t,ysol]=ode45(f,[-10 10],[2 0]);
% [t,ysol]=ode45(f,[-10 10],[0 2]);



figure;
plot(t,ysol(:,1),'r-O',t,ysol(:,2),'b-*');
title('Vander pol equation');
xlabel('x');
ylabel('y solution');
legend('y1','y2');