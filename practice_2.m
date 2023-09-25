clc;
clear;

% A=[-4 14 0;-5 13 0;-1 0 2];
% disp(A)
% tol=10^-5;
% x0=[1,1,1]';
% 
% for i=1:100
%     x1=A*x0;
%     [val,pos]=max(abs(x1));
%     eigen_val=x1(pos);
%     x2=x1/eigen_val;
%     if abs(max(x2-x1))<=tol
%         break;
%     else
%         x0=x2;
%         eigen_vec=x2;
%     end
% end
% 
% disp(eigen_val);
% disp(eigen_vec);


% syms t y;
% sol=dsolve('Dy=y*t^3 - 1.5*y','y(0)=1');
% fplot(sol);
% axis([0 2 -.5 4]);
% hold on;
% 
% f=@(t,y) y*t^3 - 1.5*y;
% y0=1;
% [tsol,ysol]=ode23(f,[0,2],y0);
% plot(tsol,ysol,'ro');


f=@(t,y) [1.2*y(1)-.6*y(1)*y(2); -.8*y(2)+.3*y(1)*y(2)];

[tsol,ysol]=ode23(f,[-10,10],[2,1]);

mu=3;
f=@(t,y) [y(2); mu*(1-y(1)^2)*y(2)-y(1)];
[tsol,ysol]=ode45(f,[-10,10],[2,0]);















































