                %%% a3Q3

clc;
clear;
A=[-4 14 0;-5 13 0;-1 0 2];

n=size(A,1);
tol=0.00001;
x=[1 1 1]';

for k=1:100
    x1=A*x;
    [val,pos]=max(abs(x1));
    eigen_val=x1(pos); 
    x2=x1/eigen_val;
    err=max(abs(x2-x)); %difference between initial guess and x2
    if err<=tol
        break
    else
        x=x2;
        eigen_vec=x;
    end
end

disp('Eigenvalue:');
disp(eigen_val);
disp('EigenVector:');
disp(eigen_vec);
fprintf('Iteration needed %d\n',k)