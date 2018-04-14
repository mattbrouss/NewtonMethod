%User Inputs. S_0 and mu control the number of outer iterations, alpha and
%beta control the linesearch function, while epsilon and maxiter control
%stopping conditions.
n = 5; % number of vars
m = 1; % number of equality constraints
condnum = 2; % (approximate) condition number for matrix Q
s_0=1;
mu=2;
epsilon=10^(-5);
maxiter=2000;
alpha=.4;
beta=.8;






Q = sprandsym(n,1.0,1/condnum,1);
c = randn(n,1);
xfeas = rand(n,1);
A = randn(m,n);
% This is to ensure that {x | Ax=b, x>=0} is nonempty
%   and there is a strictly feasible point xfeas
b = A*xfeas;

%Because Quadprog requires the optimization toolbox, we will use Yalmip and
%SDP3 as our test optimizer.
% z=sdpvar(n,1);
% constraints=[A*z-b==0,z>0];
% optimize(constraints,1/2*z'*Q*z+c'*z);

fid=fopen('OptimizationOutput.txt','w');

F=@(x) 1/2*x'*Q*x+c'*x;
GradF=@(x) Q*x+c;
HessF=@(x) Q;




[Inneriter,Solutionfval,Solutioniter,Solutionx]=Barrier(F,GradF,HessF,xfeas,A,s_0,mu,epsilon,alpha,beta,maxiter,Q,c,fid);
% zval=value(z);
% zopt=F(zval);
% Inneriter
% Solutioniter
fprintf(fid,'Optimal x=\n',Solutionx);
fprintf(fid,'%f\n',Solutionx);
fprintf(fid,'Optimal F(x)=%f\n',Solutionfval);
fprintf(fid,'Number of outer iterations=%f\n',Solutioniter);
fprintf(fid,'Total Newton Steps=%f \n',Inneriter);
fprintf(fid,'Compare to \n');
fprintf(fid,'Yalmip Optimal x=\n');
fprintf(fid,'%f\n',zval);
fprintf(fid,'Yalmip Optimal fval=%f\n',zopt);
close(fid);