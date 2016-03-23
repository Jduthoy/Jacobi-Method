function [ x ] = Jacobi1D(A,b,x0,M,Tol)
% 1D Jacobi where x is the numerical approximation for the linear system
% Ax=b, A is an nxn matrix, b nx1 vector, x0 is an initial solution, M is
% the max number of iteration, and Tol is the machine precision
[n,m]=size(A); % n=#rows m=#cols
x=zeros(n,1);
xx=A\b % answer using matlab command


for k=1:M
    r=norm(x0-xx);
       
    if r<Tol
        break;
    end
    
        x(1)=[b(1)-A(1,2:m)*x0(2:m)]/A(1,1); %Sets value for first approximation
        for i=1:n-1
            x(i)=[b(i)-(A(i,1:i-1)*x0(1:i-1)+A(i,i+1:m)*x0(i+1:m))]/A(i,i); 
        end
        x(n)=[b(n)-A(n,1:n-1)*x0(1:n-1)]/A(n,n); %Sets value for last
        x0=x;
   
 
end
r
k
end

