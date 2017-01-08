function [x,Elist,Rlist] = perform_l1_reg_fb(y,Phi,lambda, options)

% perform_l1_reg_fb - solve l1 norm regularization with FB
%
%   [x,Elist] = perform_nucl_reg_fb(y,Phi,lambda, options);
%
%   Copyright (c) 2015 Gabriel Peyre


[P,N] = size(Phi);

options.null = 0;
tau = getoptions(options, 'tau', .5 * 2/norm(Phi)^2 );
niter = getoptions(options, 'niter', 50); 
repport = getoptions(options, 'repport', @(x)0);

Thresh = @(x,t)max(1-t./max(abs(x),1e-15),0).*x;

x = zeros(N,1);
PhiS = @(z)Phi'*z;
Elist = [];
for i=1:niter
    % gradient step
    x = x - tau*PhiS( Phi*x-y );
    x = Thresh(x,lambda*tau);
    % repporting
    Elist(i) = 1/2*norm(Phi*x-y)^2+lambda*sum(abs(x(:)));
    Rlist(i) = repport(x);
end

end

