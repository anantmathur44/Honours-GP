function [l,sigma] = likelihood(x,z,par,mu,n)

if par<=0
    l = -Inf;
else
    cov = generate_cov_row(x,par);
    cov = toeplitz(cov);
    L=chol(cov,'lower'); zi=L\(z-mu);
    %sigma = (zi'*zi)/n;
    sigma = 1;
    %l = -sum(log(diag(L)))-0.5*n*log(sigma);
    l = -sum( log(diag(L)) )-0.5*(zi'*zi);
end

end