Q=eye(5);
A=[1 2 3 4 5];
c=[0; 0; 0; 0; 0];
xfeas=rand(5,1);
b=A*xfeas;
s_0=1;
mu=2;
epsilon=10^(-5);
maxiter=2000;
alpha=.4;
beta=.8;
fid=fopen('TestOutput.txt','w');

F=@(x) 1/2*x'*Q*x+c'*x;
GradF=@(x) Q*x+c;
HessF=@(x) Q;

%%This runs Barrier.m and below on the above problem, solving the problem
%%via barrier method
[Inneriter,Solutionfval,Solutioniter,Solutionx]=Barrier(F,GradF,HessF,xfeas,A,s_0,mu,epsilon,alpha,beta,maxiter,Q,c,fid);
Inneriter
Solutioniter %total inner and outer iterations, respectively
Solutionx %optimal x value

%%This runs only ConsNewton.m and below, running a constrained Newton's
%%Method on one approximation function
[Outputstatus,Outputiter,Outputfval,Outputx]=ConsNewton(F,GradF,HessF,xfeas,A,epsilon,alpha,beta,maxiter,s_0,fid);
Outputiter %total Newton steps
Outputx %optimal x value

%%This runs only Linesearch.m.
F2=@(x) s_0*F(x)-sum(log(x));
GradF2=@(x) s_0*GradF(x)-1./x;
HessF2=@(x) s_0*HessF(x)+diag(1./x.^2);
H=HessF2(x);
    G=GradF2(x);
M=[H A';A zeros(min(size(A)))];
    B=[-G; zeros(min(size(A)),1)];
    v=M\B;
    d_nt=v(1:length(x),1);
    t=Linesearch(F2,x,d_nt,G'*d_nt,alpha,beta)