function [ ite,r ] = Jacobi2D( F,M,Tol,N )
% Jacobi  with natural ordering on 2D Poisson equation where F is the
% discretization of the function f(x,y) in -(Uxx+Uyy)=f(x,y), M is the max
% number of iteration, Tol is the machine precision, N is the number of
% interior grid points with step size h=1/(N+1), and ite is the number of
% iterations the function needs for convergence. 
%syms x y
%f=-(diff(F,x,2)+diff(F,y,2));

h=1/(N+1);
U=zeros(N+2,N+2);
V1=zeros(N+2,N+2);
V2=zeros(N+2,N+2);
V3=zeros(N+2,N+2);
r=zeros(3,1);
ite=zeros(3,1);

 %{   
for i=1:N+2
   for j=1:N+2
       a=i-1;
       b=j-1;
       V(i,j)=subs(f,{x,y},{a*h,b*h});
   end
end
%}



for k=1:M
    
for i=2:N+1
    for j=2:N+1
        V1(i,j)=[U(i-1,j)+U(i+1,j)+U(i,j-1)+U(i,j+1)+(h^2)*F(i,j)]/4;
    end
end
  r(1)= norm(U-V1);
    U = V1; 
    if r(1) < Tol
       break
    end
    ite(1)=k;

end


U=zeros(N+2,N+2);
for k=1:M
    for i=2:N+1
        for j=2:N+1
            V2(i,j)=[U(i-1,j)+U(i+1,j)+U(i,j-1)+U(i,j+1)+(h^2)*F(i,j)]/4;
        end
    end
    r(2)=norm(U-V2);
    U=V2;
    if r(2)< Tol
       break
    end
    ite(2)=k;
end
  

U=zeros(N+2,N+2);

for k=1:M
    for i=2:N+1
        for j=2:N+1
            V3(i,j)=[U(i-1,j)+U(i+1,j)+U(i,j-1)+U(i,j+1)+(h^2)*F(i,j)]/4;
        end
    end
    r(3)=norm(U-V3,inf);
    U=V3; 
    if r(3)<Tol
       break
    end
    ite(3)=k;

end

