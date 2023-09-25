                %%% a3Q2

%%% Gaussian elimination mathod
clc;
clear;
A=[2 -3 2 5;-4 2 -6 14;2 2 4 8];
disp('The system is: ');
disp(A);
disp('solving the system using Gaussina elimination method');

n=size(A,1);

for j=1:n-1
    for k=2:n
        if(A(j,j)==0)
            t=A(j,:);
            A(j,:)=A(k,:);
            A(k,:)=t;
        end
    end
    for i=j+1:n
        A(i,:)=A(i,:)-A(j,:)*(A(i,j)/A(j,j));
    end
end

x=zeros(1,n);
for i=n:-1:1
    c=0;
    for k=2:n
        c=c+A(i,k)*x(k);
    end
    x(i)=(A(i,n+1)-c)/A(i,i);
end
disp('The solution is: ');
disp(x');

%%% Gauss-Jordan Elimination Method
disp('solving the system using Gauss-Jordan elimination method');
for j=1:n-1
    for k=2:n
        if(A(j,j)==0)
            t=A(j,:);
            A(j,:)=A(k,:);
            A(k,:)=t;
        end
    end
    for i=j+1:n
        A(i,:)=A(i,:)-A(j,:)*(A(i,j)/A(j,j));
    end
end

for j=n:-1:2
    for i=j-1:-1:1
        A(i,:)=A(i,:)-A(j,:)*(A(i,j)/A(j,j));
    end
end

for i=1:n
    A(i,:)=A(i,:)/A(i,i);
end

x=zeros(1,n);
for i=1:n
    x(i)=A(i,n+1);
end
disp('the solution is:');
disp(x');