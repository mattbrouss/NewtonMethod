function [Inneriter,Solutionfval,Solutioniter,Solutionx]=Barrier(F,GradF,HessF,x_0,A,s_0,mu,epsilon,alpha,beta,maxiter,Q,c,fid)
%F is the function
%GradF is the gradient thereof
%HessF is the hessian
%x0 is an initial guess (must be feasible)
%A and b the equality constraint matrix. We need to pass it down to
%ConsNewton.m
%s_0 dictates the speed/accuracy tradeoff for the initial iteration. High
%s_0 gives a slower but more accurate approximation.
%mu is the rate of centering step increase. High mu decreases the number of
%outer iterations but increases the number of inner iterations.
%epsilon is the tolerence
%alpha and beta are parameters to pass on to linesearch
%maxiter is the number of iterations you're willing to sit through before
    %giving up
%fid is the file id of the file we'll write our data to.
    
%Inneriter totals the number of Newton iterations used.  
Inneriter=0;    
x=x_0;
s=s_0;
m=(min(size(A)));
%Each iteration finds the optimal x for an approximation function, using
%that as the starting x for the next approximation.
for j=1:maxiter
    [Outputstatus,Outputiter,Outputfval,Outputx]=ConsNewton(F,GradF,HessF,x,A,epsilon,alpha,beta,maxiter,s,fid);
    x=Outputx;
    Inneriter=Inneriter+Outputiter;

    fval=F(x);
    if (m/s < epsilon)
        fprintf(fid,'Finished.\n');
    break
    
    else
        s=s*mu;
        fprintf(fid,'Outer iter %3d: Outer step=%5.2f, fval = %10.5e \n', j, s, fval);
    end
end
%Outputs
%Solutionfval gives the approximated optimal function value of the original
%F. Solutioniter gives the final number of outer iterations. Solutionx
%gives the approximated optimal x for the original F. Approximation
%accuracy depends on epsilon.
Solutionfval=1/2*x'*Q*x+c'*x;
Solutioniter=j;
Solutionx=x;
end


    
