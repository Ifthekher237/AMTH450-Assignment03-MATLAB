%%% a3Q4
%%% ode solver
clc;
clear;
syms t y(t);
od=diff(y,t)==t^3*y-1.5*y;
con=y(0)==1;
ysol(t)=dsolve(od,con);

figure
f=ezplot(ysol,[0,2]);
set(f,'linewidth',2,'color','b','linestyle','-.')
grid on
hold on
f=@(t,y) (t^3*y-1.5*y);

y0=1;

[t,y]=ode23(f,[0 2],y0);


plot(t,y,'rO-');
hold on;
f=@(t,y) (t^3*y-1.5*y);

[t,y]=ode45(f,[0 2],1);


plot(t,y,'k--*');
hold off;

%%% using euler's method with h=0.5
clear all;
f=@(t,y) (t^3*y-1.5*y);
h=0.5;
t=0:h:2;
n=length(t);
y(1)=1;
for i=1:n-1
    y(i+1)=y(i)+h*f(t(i),y(i));
end


figure
plot(t,y,'g-O');
hold on;

syms t y(t)
od=diff(y,t)==t^3*y-1.5*y;
con=y(0)==1;
ysol(t)=dsolve(od,con);

f=ezplot(ysol,[0,2]);
set(f,'linewidth',2,'color','b','linestyle','-.');
grid on;

%%% using Euler's method with h=0.25
clear;
f=@(t,y) (t^3*y-1.5*y);
h=0.25;
t=0:h:2;
n=length(t);
y(1)=1;
for i=1:n-1
    y(i+1)=y(i)+h*f(t(i),y(i));
end


figure;
plot(t,y,'g-O');
hold on;
syms t y(t);
od=diff(y,t)==t^3*y-1.5*y;
con=y(0)==1;
ysol(t)=dsolve(od,con);

f=ezplot(ysol,[0,2]);
set(f,'linewidth',2,'color','b','linestyle','-.')
grid on;

%%% Using the midpoint method with h=0.5
clear;
f=@(t,y) (t^3*y-1.5*y);
h=0.5;
t=0:h:2;
n=length(t);
y(1)=1;
for i=1:n-1
    k1=f(t(i),y(i));
    k2=f(t(i)+h/2,y(i)+h/2*k1);
    y(i+1)=y(i)+h*k2;
end


figure;
plot(t,y,'r-*')
hold on;
syms t y(t);
od=diff(y,t)==t^3*y-1.5*y;
con=y(0)==1;
ysol(t)=dsolve(od,con);

f=ezplot(ysol,[0,2]);
set(f,'linewidth',2,'color','b','linestyle','-.');
grid on;

%%% Using Heun's method with h=0.5
clear
f=@(t,y) (t^3*y-1.5*y);
h=0.5;
t=0:h:2;
n=length(t);
y(1)=1;
for i=1:n-1
    k1=f(t(i),y(i));
    k2=f(t(i)+h,y(i)+h*k1);
    y(i+1)=y(i)+h/2*(k1+k2);
end


figure;
plot(t,y,'r-*','Linewidth',2);
hold on;
syms t y(t);
od=diff(y,t)==t^3*y-1.5*y;
con=y(0)==1;
ysol(t)=dsolve(od,con);

f=ezplot(ysol,[0,2]);
set(f,'linewidth',2,'color','b','linestyle','-.')
grid on;

%%% using the fourth -order RK method with h=0.5
clear;
f=@(t,y) (t^3*y-1.5*y);
h=0.5;
t=0:h:2;
n=length(t);
y(1)=1;
for i=1:n-1
    k1=f(t(i),y(i));
    k2=f(t(i)+h/2,y(i)+h/2*k1);
    k3=f(t(i)+h/2,y(i)+h/2*k2);
    k4=f(t(i)+h,y(i)+h*k3);
    y(i+1)=y(i)+h/6*(k1+2*k2+2*k3+k4);
end


figure;
plot(t,y,'r-O','Linewidth',2);
hold on;
syms t y(t);
od=diff(y,t)==t^3*y-1.5*y;
con=y(0)==1;
ysol(t)=dsolve(od,con);

f=ezplot(ysol,[0,2]);
set(f,'linewidth',2,'color','b','linestyle','-.');
grid on;