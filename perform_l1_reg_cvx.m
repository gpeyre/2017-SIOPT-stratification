function [x,p,eta] = perform_l1_reg_cvx(y,Phi,lambda)

% perform_l1_reg_cvx - solves inverse problems with l1 norm regularization 
%
% Primal problem (here Phi(x)=Phi*x(:), x is a square matrix)
%     min |x|_1 + 1/(2*lambda)*|Phi(x)-y|^2    if   lambda>0
%     min |x|_1  s.t.  Phi(x)=y    if   lambda==0
% Dual problem
%     min |y/lambda-p|^2 s.t. |Phi^*(p)|_inf<=1    if   lambda>0
%     min <y,p>   s.t. |Phi^*(p)|_inf<=1    if   lambda==0
% and eta=Phi^*(x).
%
%   Copyright (c) 2015 Gabriel Peyre

[P,N] = size(Phi);

cvx_solver sdpt3 % SeDuMi %
cvx_begin quiet
% cvx_precision high;
variable x(N,1);
variable q(P,1);
dual variable p;
if lambda==0
    minimize norm(x,1);
    subject to
        Phi*x==y : p;
else    
    minimize norm(x,1) + 1/(2*lambda)*norm(q, 'fro')^2;
    subject to
        Phi*x-y-q==0 : p;
end
cvx_end
eta = Phi'*p;

end