function [x,Elist,Rlist] = perform_nucl_reg_fb(y,Phi,lambda, options)

% perform_fb_nucl - solve nuclear norm regularization with FB
%
%   [x,Elist] = perform_nucl_reg_fb(y,Phi,lambda, options);
%
%   Copyright (c) 2015 Gabriel Peyre


[P,N] = size(Phi);
n = sqrt(N);

options.null = 0;
tau = getoptions(options, 'tau', .5 * 2/norm(Phi)^2 );
niter = getoptions(options, 'niter', 50); 
repport = getoptions(options, 'repport', @(x)0);

Thresh = @(x,t)max(1-t./max(abs(x),1e-15),0).*x;

x = zeros(n);
PhiS = @(z)reshape(Phi'*z,[n n]);
Elist = [];
for i=1:niter
    % gradient step
    x = x - tau*PhiS( Phi*x(:)-y );
    [U,S,V] = svd(x); S = diag(S);
    S = Thresh(S,lambda*tau);
    % threshold step
    x = U*diag(S)*V';
    % repporting
    Elist(i) = 1/2*norm(Phi*x(:)-y)^2+lambda*sum(svd(x));
    Rlist(i) = repport(x);
end

end

