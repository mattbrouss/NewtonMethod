function [Outputstatus,Outputiter,Outputfval,Outputx]=ConsNewton(F,GradF,HessF,xinit,A,epsilon,alpha,beta,maxiter,s,fid)
%F is the function, GradF its gradient, abd HessF its hessian
%xinit gives the starting point,A the constraint matrix, and epsilon is tolerance
%alpha, beta line search parameters
%maxiter is maximum iterations
%s is a parameter from Barrier.m controling the approximation functions.
%Based on s_0 and mu inputs.
%fid is the file id to which we will print our results
x=xinit;
%Because this problem has constraints, we use approximations of our
%function to speed up our final iteration. These approximations are set in
%each iteration of Newton's Method based on s.
F2=@(x) s*F(x)-sum(log(x));
GradF2=@(x) s*GradF(x)-1./x;
HessF2=@(x) s*HessF(x)+diag(1./x.^2);
%Here we go through the standard Newton's Method on the approximation
%function.
for i=1:(maxiter)
    H=HessF2(x);
    G=GradF2(x);
    %Because we have equality constraints, we must  calculate our search
    %direction differently. We solve the KKT matrix and take the proper
    %length vector as our search direction (discarding the dual).
    M=[H A';A zeros(min(size(A)))];
    B=[-G; zeros(min(size(A)),1)];
    v=M\B;
    d_nt=v(1:length(x),1);
    dec=sqrt(d_nt'*H*d_nt);
    %stop condition
    if (dec^2/2 < epsilon)
        fprintf(fid,'Stop.\n');
        break;
    else
        t=Linesearch(F2,x,d_nt,G'*d_nt,alpha,beta);
        

        x=x+t*d_nt;

        fval=F(x);
        fprintf(fid,'iter %3d: step=%5.2f, fval = %10.5e \n', i, t, fval);
    end
end
%Outputs
%Outputstatus tells us if the function converged sufficiently in the number
%of iterations we allowed.
%Outputiter is the number of inner iteration for the associated outer
%iteration. Outputfval gives the optimal function value calculated for the
%approximation function. Outputx gives the calculated optimal x. These
%return to Barrier.m
if i>=maxiter
    Outputstatus='All is lost.';
    fprintf('All is lost.\n')
    fprintf(fid,'Newton Error');
else
    Outputstatus='Converged.';
end
fval=F(x);

Outputiter=i;
Outputfval=fval;
Outputx=x;
end

    

    