n = 50; % number of vars
m = 5; % number of equality constraints
condnum = 2; % (approximate) condition number for matrix Q
s_0=1;
mu=2;
epsilon=10^(-5);
maxiter=2000;
beta=.8;






Q = sprandsym(n,1.0,1/condnum,1);
c = randn(n,1);
xfeas = rand(n,1);
A = randn(m,n);
% This is to ensure that {x | Ax=b, x>=0} is nonempty
%   and there is a strictly feasible point xfeas
b = A*xfeas;
fid=fopen('Alphatest.txt','w');
X=zeros(100,1);
Y=zeros(100,1);

for i=1:100
    alpha=.5/100*i;
    F=@(x) 1/2*x'*Q*x+c'*x;
GradF=@(x) Q*x+c;
HessF=@(x) Q;
[Inneriter,Solutionfval,Solutioniter,Solutionx]=Barrier(F,GradF,HessF,xfeas,A,s_0,mu,epsilon,alpha,beta,maxiter,Q,c,fid);
X(i)=alpha;
Y(i)=Inneriter;
end
figure
plot(X,Y)
title('Convergence via Alpha')
xlabel('Alpha')
ylabel('Total Iterations')

alpha=.2;
for i=1:100
    beta=.9999/100*i;
    F=@(x) 1/2*x'*Q*x+c'*x;
GradF=@(x) Q*x+c;
HessF=@(x) Q;
[Inneriter,Solutionfval,Solutioniter,Solutionx]=Barrier(F,GradF,HessF,xfeas,A,s_0,mu,epsilon,alpha,beta,maxiter,Q,c,fid);
X(i)=beta;
Y(i)=Inneriter;
end
figure
plot(X,Y)
title('Convergence via Beta')
xlabel('Beta')
ylabel('Total Iterations')

beta=.8;

for i=1:100
    s_0=i/5;
    F=@(x) 1/2*x'*Q*x+c'*x;
GradF=@(x) Q*x+c;
HessF=@(x) Q;
[Inneriter,Solutionfval,Solutioniter,Solutionx]=Barrier(F,GradF,HessF,xfeas,A,s_0,mu,epsilon,alpha,beta,maxiter,Q,c,fid);
X(i)=s_0;
Y(i)=Inneriter;
end
figure
plot(X,Y)
title('Convergence via S_0')
xlabel('S_0')
ylabel('Total Iterations')
s_0=1;

for i=1:100
    mu=1+10*i;
    F=@(x) 1/2*x'*Q*x+c'*x;
GradF=@(x) Q*x+c;
HessF=@(x) Q;
[Inneriter,Solutionfval,Solutioniter,Solutionx]=Barrier(F,GradF,HessF,xfeas,A,s_0,mu,epsilon,alpha,beta,maxiter,Q,c,fid);
X(i)=mu;
Y(i)=Inneriter;
end
figure
plot(X,Y)
title('Convergence via mu')
xlabel('mu')
ylabel('Total Iterations')
mu=2;


% This is to ensure that {x | Ax=b, x>=0} is nonempty
%   and there is a strictly feasible point xfeas
for i=1:100
    n=1+10*i;
    Q = sprandsym(n,1.0,1/condnum,1);
c = randn(n,1);
xfeas = rand(n,1);
A = randn(m,n);
b = A*xfeas;
F=@(x) 1/2*x'*Q*x+c'*x;
GradF=@(x) Q*x+c;
HessF=@(x) Q;
[Inneriter,Solutionfval,Solutioniter,Solutionx]=Barrier(F,GradF,HessF,xfeas,A,s_0,mu,epsilon,alpha,beta,maxiter,Q,c,fid);
X(i)=n;
Y(i)=Inneriter;
end
figure
plot(X,Y)
title('Convergence via n')
xlabel('n')
ylabel('Total Iterations')
n=50;

X=zeros(50,1);
Y=zeros(50,1);
for i=1:50
    m=i;
    Q = sprandsym(n,1.0,1/condnum,1);
c = randn(n,1);
xfeas = rand(n,1);
A = randn(m,n);
b = A*xfeas;
F=@(x) 1/2*x'*Q*x+c'*x;
GradF=@(x) Q*x+c;
HessF=@(x) Q;
[Inneriter,Solutionfval,Solutioniter,Solutionx]=Barrier(F,GradF,HessF,xfeas,A,s_0,mu,epsilon,alpha,beta,maxiter,Q,c,fid);
X(i)=m;
Y(i)=Inneriter;
end
figure
plot(X,Y)
title('Convergence via m')
xlabel('m')
ylabel('Total Iterations')
m=5;

X=zeros(10,1);
Y=zeros(10,1);
for i=1:10
    condnum=i;
    Q = sprandsym(n,1.0,1/condnum,1);
c = randn(n,1);
xfeas = rand(n,1);
A = randn(m,n);
b = A*xfeas;
F=@(x) 1/2*x'*Q*x+c'*x;
GradF=@(x) Q*x+c;
HessF=@(x) Q;
[Inneriter,Solutionfval,Solutioniter,Solutionx]=Barrier(F,GradF,HessF,xfeas,A,s_0,mu,epsilon,alpha,beta,maxiter,Q,c,fid);
X(i)=condnum;
Y(i)=Inneriter;
end
figure
plot(X,Y)
title('Convergence via condnum')
xlabel('condnum')
ylabel('Total Iterations')



