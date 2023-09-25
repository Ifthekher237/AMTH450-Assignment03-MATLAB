                %%% a3Q1

%%% ans to Q no a
%Jacobi iterative method
disp('Jacobi iterative method')
clc;
clear;
A=[10,-1,2,0;-1,11,-1,3;2,-1,10,-1;0,3,-1,8];
b=[6;25;-11;15];
n=size(A,1);
x0=[0 0 0 0];
x1=x0;
kmax=30;
tol=0.00001;
k=1;
while k<=kmax
    for i=1:n
        s=0;
        for j=1:n
            if j~=i
                s=s+A(i,j)*x0(j);
            end
        end
    s=(b(i)-s)/A(i,i);
    x1(i)=s;
    end
    err=max(abs(x1-x0));
    if err<=tol
        break
    else
        k=k+1;
        x0=x1;
    end
end
fprintf('Iteration needed: %d\n',k)

for i=1:n
    fprintf('Value of x%d is %11.5f\n',i,x1(i))
end

disp(' ')
%%% ans to Q no b
%Gauss-Seidel iterative method
disp('Gausss-Seidel iterative method')
A=[10,-1,2,0;-1,11,-1,3;2,-1,10,-1;0,3,-1,8];
b=[6;25;-11;15];
n=size(A,1);
x0=[0 0 0 0];
x1=x0;
kmax=30;
tol=0.00001;
k=1;
while k<=kmax
    for i=1:n
        s=0;
        for j=1:n
            if j~=i
                s=s+A(i,j)*x1(j);
            end
        end
        s=(b(i)-s)/A(i,i);
        x1(i)=s;
    end
    err=max(abs(x1-x0));
    if err<=tol
        break
    else
        k=k+1;
        x0=x1;
    end
end
fprintf('Iteration needed: %d\n',k)

for i=1:n
    fprintf('Value of x%d is %11.5f\n',i,x1(i))
end

disp(' ')
%%% ans to q no c
%SOR iterative method with w=1.1
disp('SOR iterative method with w=1.1')
A=[10,-1,2,0;-1,11,-1,3;2,-1,10,-1;0,3,-1,8];
b=[6;25;-11;15];
n=size(A,1);
x0=[0 0 0 0];
x1=x0;
w=1.1;
kmax=30;
tol=0.00001;
k=1;
while k<=kmax
    for i=1:n
        s=0;
        for j=1:n
            if j~=i
                s=s+A(i,j)*x1(j);
            end
        end
        s=(w/A(i,i))*(b(i)-s);
        x1(i)=(1-w)*x0(i)+s;
    end
    err=max(abs(x1-x0));
    if err<=tol
        break
    else
        k=k+1;
        x0=x1;
    end
end
fprintf('Iteration needed: %d\n',k)

for i=1:n
    fprintf('Value of x%d is %11.5f\n',i,x1(i))
end