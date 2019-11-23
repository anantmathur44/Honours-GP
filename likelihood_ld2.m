function [l,g] = likelihood_ld2(x,z,par,mu,n)

if (par(1)<=0) || (par(2) <=0)
    l = -Inf;
    if nargout > 1 % gradient required
        g = -mc_score_var(x,z,n-1,[par(1), par(2)],10);     
    end
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
    k = -sum( log(diagon) )-0.5*((U-v*mu)'*(U-v*mu));
    l = -k;
    if nargout > 1 % gradient required
        g = -mc_score_var(x,z,n-1,[par(1), par(2)],10);     
    end
    
end

end

