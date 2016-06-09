function [p,eta] = compute_certificate_nucl(x,Phi)

% compute_certificate_nucl - compute minimal norm certificate
%
%   [p,eta] = compute_certificate_nucl(x,Phi);
%
%   Copyright (c) 2015 Gabriel Peyre

[P,N] = size(Phi);
n = sqrt(N);


% identify support
[U,S,V] = svd(x); S = diag(S);
tol = 1e-5;
I = find(S/max(S)>=tol);
J = find(S/max(S)<tol);
U1 = U(:,J); V1 = V(:,J); 
e = U(:,I)*V(:,I)';

% min |p|^2  s.t. norm(
cvx_solver sdpt3 % SeDuMi %
cvx_begin quiet
% cvx_precision high;
variable p(P,1);
variable eta(n,n);
minimize norm(p,'fro');
subject to
    Phi'*p-eta(:)==0;
    %P_T(eta)=e
    eta-(U1*U1')*eta*(V1*V1')-e==0;
    norm(eta)<=1; 
cvx_end

end