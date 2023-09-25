%%% a3Q6

clear;
clc;
p=1;
cp=1;
k=1;
alpha=k/(p*cp);
xmin=0;
xmax=10;
N=50;
dx=(xmax-xmin)/(N-1);
x=xmin:dx:xmax;
dt=0.001;
tmax=5;
t=0:dt:tmax;
s(1:N,1)=ones(1,N)*100;
ux0=0;
ux10=ux0;
r=alpha*dt/(dx^2);
for j=2:length(t)
    u=s(1:N,j-1);
    for i=1:N
        if i==1|| i==N
            s(i,j)=ux0;
        else
            s(i,j)=u(i)+r*(u(i+1)-2*u(i)+u(i-1));

        end
    end
end


figure;
imagesc('XData',t,'YData',x,'CData',s);
xlabel('time t');
ylabel('space x');
colorbar;
title('Graphical visualization of heat conduction with time and space');