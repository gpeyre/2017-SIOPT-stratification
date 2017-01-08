function [p,eta] = compute_certificate_l1(x,Phi)

% compute_certificate_l1 - compute minimal norm certificate
%
%   [p,eta] = compute_certificate_l1(x,Phi);
%
%   Copyright (c) 2015 Gabriel Peyre

[P,N] = size(Phi);

% identify support
tol = 1e-5;
I = find(abs(x)>=tol);

% min |p|^2  s.t. norm(
cvx_solver sdpt3 % SeDuMi %
cvx_begin quiet
% cvx_precision high;
variable p(P,1);
variable eta(N,1);
minimize norm(p,'fro');
subject to
    Phi'*p-eta==0;
    %P_T(eta)=e
    eta(I)==sign(x(I));
    norm(eta, Inf)<=1; 
cvx_end

end