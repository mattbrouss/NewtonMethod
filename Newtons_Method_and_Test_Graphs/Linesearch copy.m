function t = Linesearch(F, x, d, df, alpha, beta, tinit)
% F     -- Function evaluation
% x     -- current x
% d     -- search direction
% df    -- directional derivative in d at t=0, i.e.,
%         it equals gradF(x)'*d
% alpha, beta -- parameters in line search
% tinit -- initial step size
% t     --  Final step size

%Initialize, check parameter bounds
if (nargin < 7)
    tinit = 1;
end

if (alpha<=0 || alpha >= 1 || ...
        beta <= 0 || beta >= 1)
    fprintf(fid,'Error: alpha and beta must be in proper range.\n');
    t = 0;
    return;
end
t = tinit;
% This section iterates t=beta*t until x+t*d is element wise positive
for i=1:length(x)
    while (x(i)+t*d(i)<0)
      t = t*beta;
    end
end
%Using t sufficiently small (as found above) we continue with standard line
%search
while (F(x+t*d) > (F(x) + alpha*t*df))
    t = t*beta;
end
%Only output is t, the final step length. Returns to ConsNewton.m
