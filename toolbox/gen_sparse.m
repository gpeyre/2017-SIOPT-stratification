function x = gen_sparse(N,s)
    
x = zeros(N,1);
I = randperm(N); I = I(1:s); 
x(I) = sign(randn(s,1));

end