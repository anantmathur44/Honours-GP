function [l,sigma] = likelihood_ld(x,z,par,mu,n)

if (par(1)<=0) || (par(2) <=0)
    l = -Inf;
else
    gam = generate_cov_row(x,par);
    %L=chol(cov,'lower'); zi=L\z;    
    [L,diagon]=LGL(gam); diagon=sqrt(diagon');
    U=(L*z)./diagon;
    v = sum(L,2)./diagon;
    %b = (v'*U)/(v'*v);
    %sigma = ((U-v*mu)'*(U-v*mu))/n;
    sigma = 1;
    %l = -sum(log(diagon))-0.5*n*log(sigma);
    l = -sum( log(diagon) )-0.5*((U-v*mu)'*(U-v*mu));
end

end

