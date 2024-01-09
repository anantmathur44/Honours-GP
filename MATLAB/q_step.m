function q = q_step(x,Y,Z,par,m,mu,n)

%generating circ row
col = generate_cov_row(x,par);
circ = [col; col(end-1:-1:2)];
prods = [];
global mu_q;
global sigma_sq_q;
if (par(1)<=0)||(par(2) <=0)
    q = -Inf;
    return
else
    %eigenvalues
    lambda = ifft(circ);
    %iterate through samples
%     mu_q = mean(mean(Z));

    
    for i = Z
        comp_z = [Y;i];
        %comp_z'*inv(toeplitz(circ))*comp_z
        quad = sqrt(lambda).^-1.*ifft(comp_z); 
        prods = [prods (norm(quad)^2)/m];
        %prods = [prods quadratic_product(lambda,comp_z,1)/n];
    end
    
%     sigma_sq_q = (sum(prods))/(2*n);
%     q = -n*log(sigma_sq_q) - 0.5*sum(log(lambda));


      q = -0.5*(sum(prods)+sum(log(lambda*2*n)));
  
end

end