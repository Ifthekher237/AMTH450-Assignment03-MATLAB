 clc;
 clear;
 
%  A=transpose(reshape([1:1:40],[8,5]));
% disp(A); 
% 
% sub_A=A([1,3,5],[1,2,4,8]);
% disp(sub_A);
% 
% sub_A2=unique([A(5,:) transpose(A(:,4)) transpose(A(:,6))],'stable');
% disp(sub_A2);

% a=.75;b=11.3;
% x=[2,5,1,9];y=[.2,1.1,1.8,2];z=[-3,2,5,4];
% 
% A1=((power(x,1.1).*power(y,-2).*power(z,5))./(power(a+b,b/3))) + a.*(((z./x)+(y./2))/(power(z,a)));
% disp(A1)

% syms x1 x2 x3 x4;
% 
% [A,B]=equationsToMatrix(2*x1+x2+x3-x4==12,x1+5*x2-5*x3+6*x4==35,-7*x1+3*x2-7*x3-5*x4==7,x1-5*x2+2*x3+7*x4==21);
% X=linsolve(A,B);
% disp(X);
% [x1,x2,x3,x4]=solve(2*x1+x2+x3-x4==12,x1+5*x2-5*x3+6*x4==35,-7*x1+3*x2-7*x3-5*x4==7,x1-5*x2+2*x3+7*x4==21,[x1,x2,x3,x4]);
% disp(x1);
% disp(x2);
% disp(x3);
% disp(x4);

% syms x;
% 
% x=0:.1:2*pi;
% y1=power(sin(x),2);
% y2=power(cos(x),2);
% y3=cos(2.*x);
% 
% figure(1)
% plot(x,y1);
% hold on;
% plot(x,y2);
% plot(x,y3);
% axis([0 2*pi, -2 2]);
% title('ploting trigonometric functions')
% legend('sin^2(x)','cos^2(x)','cos(2x)','location','southeast');
% line([0,2*pi],[0,0],'color','k','linewidth',1.2);
% line([0,0],[-2 2],'color','k','linewidth',1.2);
% grid on;
% 
% figure(2);
% subplot(3,1,1);
%     plot(x,y1);
%     hold on;
% subplot(3,1,2);
%     plot(x,y2);
% subplot(3,1,3);
%     plot(x,y3);
%     hold off;

% syms x y z f;
% f=.56.*cos(x.*y);
% [x,y]=meshgrid(0:.1:5, 0:.1:5);
% z=.56.*cos(x.*y);
% surf(x,y,z);
% contour(x,y,z);
% fcontour(f);

% function practice(a,b)
%     syms x y f;
%     f=x^2 + y^2 - 2*x*y + 4;
%     grad(x,y)=gradient(f);
%     disp(grad(a,b));
%     
% end

% syms x;
% sol=vpa(solve(x^7 - 8*x^5 + 7*x^4 + 5*x^3 - 8*x + 9==0,x),6);
% disp(sol);
% 
% odsolve=dsolve('D2x + 10*Dx + 5*x==11','x(0)==1','Dx(0)=-1');
% disp(odsolve);
% 
% F=x^5 - 8*x^4 + 5*x^3 - 7*x^2 + 11*x - 9;
% disp(diff(F,x,1));
% disp(diff(F,x,2));
% 
% f=1/(.8*x^2 + .5*x + 2);
% disp(int(f,x,0,5));

% year=1930:10:2020;
% pop=[249 277 316 350 431 539 689 833 1014 1203];
% 
% [p,s,mu]=polyfit(year,pop,2);
% new_pol=polyval(p,year,s,mu);
% 
% plot(year,pop,'ko');
% hold on;
% plot(year,new_pol);
% hold off;
% 
% val1=interp1(year,pop,1995,'linear');
% val2=interp1(year,pop,1995,'spline');
% disp(val1);
% disp(val2);

% syms x y;
% u=1+y^2;
% v=x-1;
% 
% y1=int(v,x);
% y2=int(u,y);
% s_line=simplify(y1-y2);
% disp(s_line);
% 
% [x,y]=meshgrid(-4:.5:4 -4:.5:4);
% streamslice(x,y,subs(u),subs(v));


% syms x y z;
% v=[(x^2/2 - x^3/3) x*(x-1)*(y+1) 0];
% 
% %motion is possible, if divergence is zero
% div=divergence([v(1,1) v(1,2)],[x,y]);
% if div==0;
%     disp('motion is possible');
% else
%     disp('motion is not possible');
% end
% 
% %motion is irrotational, if curl is zero
% crl=simplify(curl([v(1,1) v(1,2) 0],[x,y,z]));
% if crl==0
%     disp('motion is irrotationa');
% else
%     disp('rotational')
% end
% 
% %finding the stagnation points
% solx=unique(solve(v(1,1)==0,x));
% soly=unique(solve(v(1,2)==0,y));
% [x,y]=meshgrid(solx,soly);
% 
% X=transpose([x ;y]);
% disp(X);

% syms x y z a;
% 
% u=a*(x^2 - y^2);
% v=-2*a*x*y;
% 
% %stream function exists, if divergence is zero
% div=simplify(divergence([u v],[x y]));
% disp(div);
% 
% %potential function exists, if it's irrotational
% crl=simplify(curl([u v 0],[x y z]));
% disp(crl);
% 
% pot_func=potential([u v 0],[x y z]);
% disp(pot_func);
% 
% strm_func=potential([-v u 0],[x y z]);
% disp(strm_func);
% 
% f1=subs(strm_func,a,3);
% f2=subs(pot_func,a,3);
% val=-4:2:4;
% hold on;
% fcontour(f1,'LevelList',val);
% fcontour(f2,'LevelList',val);
% hold off;


% syms x y a;
% phi=(a*x^3)/3 - (a*x*y^2) -2;
% f1=subs(phi,a,3);
% val=-4:2:4;
% fcontour(f1,'LevelList',val);
% hold on;
% u=diff(phi,x);
% v=diff(phi,y);
% 
% strm_func=potential([-v u],[x,y]);
% disp(strm_func);
% f2=subs(strm_func,a,3);
% fcontour(f2,'k-','LevelList',val);
% [x,y]=meshgrid(0:1:5, 0:1:5);
% u1=subs(u,a,3);
% v1=subs(v,a,3);
% streamslice(x,y,subs(u1),subs(v1));
% hold off;

% a=.8;
% b=6;
% z2=a*.61;
% g=9.81;
% 
% i=1;
% for z1=5:15
%     Qb(i)=z2*power((2*g*z1)/(1+(z2/z1)),1/2);
%     z(i)=z1;
%     i=i+1;
% end
% plot(z,Qb);

% for n=1:5
%     R(n,1)=n;
%     Mn=power(2,n)-1;
%     R(n,2)=Mn;
%     Fn=power(2,power(2,n))+1;
%     R(n,4)=Fn;
%     
%     if isprime(Mn)
% %         prime_Mn=Mn;
%         R(n,3)=Mn;
%     else
% %         prime_Mn=0;
%         R(n,3)=0;
%     end
%     
%     if isprime(Fn)
% %         prime_Fn=Fn;
%         R(n,5)=Fn;
%     else
% %         prime_Fn=0;
%         R(n,5)=0;
%     end 
% %     fprintf('%i %i %i %i %i\n\n',n,Mn,prime_Mn,Fn,prime_Fn);
% %     R(n,:)=[n Mn prime_Mn Fn prime_Fn];
% %     fprintf('a=%d b=%d c=%d d=%d e=%d\n',n,Mn,prime_Mn,Fn,prime_Mn);
% end
% variables={'n','Mn','prime_Mn','Fn','prime_Fn'};
% T=array2table(R,'VariableNames',variables);
% disp(T);


% n=6120;
% factors=factor(n);
% uni=unique(factors);
% disp(uni);
% k=length(uni);
% for i=1:k
%     p(i)=uni(i);
%     alpha(i)=length(find(factors==p(i)));
% end

% syms x f;
% f(x)=exp(x) - 2 -cos(exp(2)-2);
% tol=.000005;
% a=.5;
% b=1.5;
% 
% if f(a)*f(b)>0
%     disp('inappropriate initial guess');
% end
% 
% if f(a)>0
%     [a,b]=deal(b,a);
% end
% 
% for i=1:1000
%     root=(a+b)/2;
%     fprintf('i=%d,a=%f,b=%f,c=%f\n',i,a,b,root);
%     if abs(f(root))<=tol
%         break;
%     else
%         if f(a)*f(root)>0
%             a=root;
%         else
%             b=root;
%         end
%     end
%     
% end


% syms x;
% tol=.000005;
% f=@(x) x.^3 - x + exp(x);
% g{1}=@(x) x.^3 + exp(x);
% g{2}=@(x) (x-exp(x)).^(1/3);
% g{3}=@(x) log(x-x.^3);
% 
% 
% for i=1:3
%     clear x;
%     syms x;
%     df(x)=diff(g{i}(x),x);
%     x=linspace(-3.1,3,100);
%     if max(abs(df(x)))>=1
%         continue;
%     end
%     disp(g{i});
%     x0=-1;
%     for k=1:1000
%         x1=g{i}(x0);
%         if isreal(x1)==0
%             x1=-abs(x1);
%         end
%         
%         fprintf('i=%d, x0=%f, x1=%f\n',k,x0,x1);
%         if abs(x1-x0)<=tol
%             break;
%         else
%             x0=x1;
%         end
%     end
% end


% syms x;
% f=@(x) 54*x^6 + 45*x^5 - 102*x^4 - 69*x^3 + 35*x^2 + 16*x-4;
% % x=-10:.01:10;
% % plot(x,f(x));
% % hold on;
% % axis([-5 5 -25 8]);
% % hold off
% 
% tol=.000005;
% solutions=unique(vpasolve(f(x)==0,x));
% disp(solutions);
% fplot(f);
% hold on;
% plot(solutions,0,'ko');
% axis([-5 5 -25 8])
% line([-5 5],[0 0],'color','k','linewidth',1.2)
% line([0 0],[-25 8],'color','k','linewidth',1.2)
% hold off;
% 
% x0=[-2,-1,0,0.45,1];
% x1=[-3,0,.25,1,2];
% % x0 = [-3 -3 0 0.45 3];
% % x1 = [-2 -1 0.25 1 2];
% 
% for i=1:5
%     disp(solutions(i));
%         p0=x0(i);
%         p1=x1(i);
%     for k=1:200
%         p=(p1*f(p0) - p0*f(p1))/(f(p0) - f(p1));
%         fprintf('i=%d, p0=%f, p1=%f, p2=%f\n',k,p0,p1,p);
%         if abs(p-solutions(i))<tol
%             break;
%         end
%         
%         p0=p1;
%         p1=p;
%     end
% end

% syms x;
% f=@(x) 14.*x.*exp(x-2) - 12.*exp(x-2) - 7.*x.^3 + 20.*x.^2 - 26.*x + 12;
% sol=[vpasolve(f(x)==0,x,1) vpasolve(f(x)==0,x,3)];
% disp(sol);
% fplot(f);
% hold on;
% axis([-2 4 -2 10]);
% line([-2 4],[0 0],'color','k','linewidth',1.2);
% line([0 0],[-2 10],'color','k', 'linewidth',1.2);
% plot(sol,0,'ko');
% grid on;
% hold off;
% 
% x0=[1,1.5];
% df(x)=diff(f(x),x,1);
% 
% for i=1:2
%     disp(sol(i));
%     p0=x0(i);
%     for j=1:200
%         p=double(p0-(f(p0)/df(p0)));
%         fprintf('i=%d, p0=%f, p=%f\n',j,p0,p);
%         if abs(p-p0)<=.000005
%             break;
%         else
%             p0=p;
%         end
%     end
% end


% A=[10 -1 2 0;-1 11 -1 3; 2 -1 10 -1;0 3 -1 8];
% B=[6;25;-11;15];
% tol=.00001;
% n=size(A,1);
% 
% x0=[0 0 0 0];
% x1=x0;
% w=1.1;
% for k=1:20
%     for i=1:n
%         s=0;
%         for j=1:n
%             if j~=i
%                 s=s+A(i,j)*x1(j);
%             end
%         end
%         s=(w*(B(i)-s))/A(i,i);
%         x1(i)=s+(1-w)*x0(i);
%     end
%     if max(abs(x1 - x0))<=tol
%         break
%     else
%         x0=x1;
%     end
% end
% disp(k);
% disp(x1);






% A=[2 -3 2 5;-4 2 -6 14;2 2 4 8];
% disp(A);
% n=size(A,1);
% 
% for i=1:n-1
%     for j=i+1:n
%         A(j,:)=A(j,:) - (A(j,i)/A(i,i))*A(i,:);
%     end
% end
% 
% for i=n:-1:1
%     c=0;
%     for j=i+1:n
%         c=c+A(i,j)*x(j);
%     end
%     x(i)=(A(i,n+1)-c)/A(i,i);
% end
% disp(x);
% 
% %for gauss-jordan
% 
% for i=n:-1:2
%     for j=i-1:-1:1
%         A(j,:)=A(j,:) - (A(j,i)/A(i,i))*A(i,:);
%     end
% end
% disp(A);
% for i=1:n
%     A(i,:)=A(i,:)/A(i,i);
% end
% disp(A);
% 
% for i=1:n
%     x(i)=A(i,n+1);
% end
% disp(x);

% tol=.00001;
% A=[-4 14 0;-5 13 0;-1 0 2];
% x0=[1 1 1]';
% 
% for i=1:200
%     x1=A*x0;
%     [val,pos]=max(abs(x1));
%     eigen_val=x1(pos);
%     x2=x1/eigen_val;
%     err=max(abs(x2 - x0));
%     if err<=tol
%         break;
%     else
%         x0=x2;
%         eigen_vec=x0;
%     end
% end
% 
% disp(eigen_vec);
% disp(eigen_val);
% disp(i);


% syms t y(t);
% od=diff(y,t)==y*t^3 - 1.5*y;
% con=y(0)==1;
% sol=dsolve(od,con);
% disp(sol);
% t=0:.001:2;
% plot(t,subs(sol));
% hold on;
% 
% % ode23
% f=@(t,y) y.*t.^3 - 1.5.*y;
% y0=1;
% [t,y]=ode23(f,[0,2],y0);
% 
% plot(t,y,'b*');
% 
% % ode 45
% f=@(t,y) y.*t.^3 - 1.5.*y;
% y0=1;
% [t,y]=ode45(f,[0,2],y0);
% 
% plot(t,y,'r*');
% legend('exact sol','ode23 results','ode45 results');
% hold off;

% eulers method with h=.5 and .25

% f=@(t,y) y.*t.^3 - 1.5.*y;
% 
% h=.5;
% t=0:h:2;
% y=1;
% 
% for i=1:length(t)-1
%     y(i+1)=y(i)+h*f(t(i),y(i));
% end
% 
% plot(t,y,'r-*');
% hold on;
% 
% %midpoint method with h=.5
% f=@(t,y) y.*t.^3 - 1.5.*y;
% t=0:.5:2;
% y=1;
% 
% for i=1:length(t)-1
%     k1=f(t(i),y(i));
%     k2=f(t(i)+h/2, y(i)+(k1*h)/2);
%     y(i+1)=y(i)+k2*h;
% end
% plot(t,y,'k--o')
% 
% %heuns method with h=.5
% f=@(t,y) y.*t.^3 - 1.5.*y;
% t=0:.5:2;
% y=1;
% 
% for i=1:length(t)-1
%     y_p(i+1)=y(i) + h*f(t(i),y(i));
%     y(i+1)=y(i) + (h/2)*(f(t(i),y(i))+f(t(i+1),y_p(i+1)));
% end
% plot(t,y,'c-+');


% Ans to Q no 5(a)
% 
% f=@(t,x) [1.2*x(1) - .6*x(1)*x(2); -.8*x(2) + .3*x(1)*x(2)];
% [t,xsol]=ode23(f,[-10,10],[2,1]);
% plot(t,xsol(:,1),'r-*',t,xsol(:,2),'b-o');
% 
% % Ans to Q no 5(b)
% mu=3;
% f=@(t,y) [y(2);mu*(1-y(1)^2)*y(2)-y(1)];
% [t,ysol]=ode45(f,[-10,10],[2,0]);
% plot(t,ysol(:,1),'r-o',t,ysol(:,2),'k-*');


% A=[10 -1 2 0;-1 11 -1 3; 2 -1 10 -1;0 3 -1 8];
% B=[6;25;-11;15];
% tol=.00001;
% n=size(A,1);
% x0=[0 0 0 0];
% x1=x0;
% w=1.1;
% for k=1:1000
%     for i=1:n
%         s=0;
%         for j=1:n
%             if j~=i
%                 s=s+A(i,j)*x1(j);
%             end
%         end
%         s=(-s+B(i))/A(i,i);
%         x1(i)=(1-w)*x0(i) + (s*w);
%     end
%     if max(abs(x1-x0))<=tol
%         break;
%     else
%         x0=x1;
%     end
% end
% disp(k)
% disp(x1)


% A=[2 -3 2 5;-4 2 -6 14; 2 2 4 8];
% n=size(A,1);
% 
% for i=1:n-1
%     for j=i+1:n
%         A(j,:)=A(j,:)-(A(j,i)/A(i,i))*A(i,:);
%     end
% end
% disp(A);
% %backward sbustitutions
% 
% % for i=n:-1:1
% %     c=0;
% %     for j=i+1:n
% %         c=c+A(i,j)*x(j);
% %     end
% %     x(i)=(A(i,n+1)-c)/A(i,i);
% % end
% % 
% % disp(x);
% 
% for i=n:-1:2
%     for j=i-1:-1:1
%         A(j,:)=A(j,:)-(A(j,i)/A(i,i))*A(i,:);
%     end
% end
% 
% for i=1:n
%     A(i,:)=A(i,:)/A(i,i);
% end
% disp(A)
% for i=1:n
%     x(i)=A(i,n+1);
% end
% disp(x);


% A=[-4 14 0;-5 13 0;-1 0 2];
% x0=[1 1 1]';
% tol=.00001;
% 
% for i=1:100
%     x1=A*x0;
%     [val,pos]=max(abs(x1));
%     eigen_val=x1(pos);
%     x2=x1/eigen_val;
%     if max(abs(x2-x0))<=tol
%         break;
%     else
%         x0=x2;
%         eigen_vec=x2;
%     end
% end
% 
% disp(i);
% disp(eigen_vec);
% disp(eigen_val);


% 
% sol=dsolve('Dy=y*t^3 - 1.5*y','y(0)=1');
% disp(sol);
% 
% figure;
% fplot(sol);
% hold on;
% axis([0 2 -1 5]);
% 
% f=@(t,y) y*t^3 - 1.5*y;
% y0=1;
% [tsol,ysol] = ode23(f,[0,2],y0);
% plot(tsol,ysol,'r-o');
% 
% 
% [tsol,ysol] = ode45(f,[0,2],y0);
% plot(tsol,ysol,'k-*');
% 
% t1=0:.5:2;
% y1=1;
% 
% for i=1:length(t1)-1
%     y1(i+1)=y1(i)+.5*f(t1(i),y1(i));
% end
% 
% t2=0:.25:2;
% y2=1;
% 
% for i=1:length(t2)-1
%     y2(i+1)=y2(i)+.25*f(t2(i),y2(i));
% end
% 
% plot(t1,y1,'b-+',t2,y2,'g-s');
% 
% t=0:.5:2;
% y=1;
% for i=1:length(t)-1
%     k1=f(t(i),y(i));
%     k2=f(t(i)+(1/2)*.5, t(i)+(1/2)*k1*.5);
%     y(i+1)=y(i)+k2*.5;
% end
% 
% plot(t,y,'y-d');
% 
% 
% 
% legend('exact','ode23','ode45','euler(h=.5)','euler(h=.25','midpoint(h=.5)');
% hold off;

% f=@(t,y) [1.2*y(1) - .6*y(1)*y(2);-.8*y(2)+.3*y(1)*y(2)];
% [tsol,ysol]=ode23(f,[-10,10],[2,1]);
% plot(tsol,ysol)
% mu=3;
% 
% f=@(t,y) [y(2);mu*(1-y(1)^2)*y(2)-y(1)];
% [tsol,ysol]=ode45(f,[-10,10],[2,0]);
% plot(tsol,ysol);

% syms x y r t;
% m=-.314;
% F=vpa((m/2*pi)*(log(x+i*y+5i)+log(x+i*y-5i)));
% psi_cart=imag(F);
% disp(psi_cart);
% psi_polar=subs(psi_cart,{x,y},{r*cos(t),r*sin(t)});
% disp(psi_polar);
% 
% grad=gradient(psi_cart,[x,y]);
% u=grad(1);
% v=-grad(2);
% disp(u);
% disp(v);
% 
% fcontour(imag(F),'b',[-10,10],'LevelList',-10:.1:10,'linewidth',1);



n=4;
if mod(n,2)==0
    N=n;
else
    N=n+1;
end

for k=1:N-1
    for i=1:N-1
        j=mod(k-i,N-1);
        
        if j==0
            j=N-1;
        end
        
        if j==i
            R(k,i)=N;
        elseif j==0
            R(k,i)=N-1;
        else
            R(k,i)=j;
        end
    end
end

disp(R);


























