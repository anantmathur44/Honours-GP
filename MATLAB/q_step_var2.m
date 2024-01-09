function q = q_step_var2(x,Y,Z,par,m,mu,n,mode,par_1)

par_previous = par_1;
par_new = par;

%generating circ row
col = generate_cov_row(x,par);
circ = [col; col(end-1:-1:2)];
circ_1 = circ;
prods = [];
prods_2 = [];
prods_3 = [];

%means = [];
global mu_q;
global sigma_sq_q;


if (par(1)<=0)||(par(2) <=0)
    q = -Inf;
    return
    
else
    %eigenvalues
    lambda = ifft(circ); 
    %iterate through samples
    %mean(mean(Z))
 
    
%   mu_q = [Y; mean(Z')'];
    mode_z = [Y;mode];
    
    for i = Z
        comp_z = [Y;i];
        %comp_z'*inv(toeplitz(circ))*comp_z
        quad = sqrt(lambda).^-1.*ifft(comp_z-mode_z); 
        prods = [prods (norm(quad)^2)/m];
        %prods = [prods quadratic_product(lambda,comp_z,1)/n];
    end
    
    
    for i = Z
        comp_z = [Y;i];
        %comp_z'*inv(toeplitz(circ))*comp_z
        quad = sqrt(lambda).^-1.*ifft(comp_z); 
        prods_3 = [prods_3 (norm(quad)^2)/m];
        %prods = [prods quadratic_product(lambda,comp_z,1)/n];
    end
    
    mu =0;
    
    quad_m = sqrt(lambda).^-1.*ifft(mode_z-mu); 
    quad_m = norm(quad_m)^2;
  
    col = generate_cov_row(x,par_1);
    circ = [col; col(end-1:-1:2)];
    lambda2 = ifft(circ); 
    
    for i = Z
        comp_z = [Y;i];
        %comp_z'*inv(toeplitz(circ))*comp_z
        quad = sqrt(lambda2).^-1.*ifft(comp_z-mode_z); 
        prods_2 = [prods_2 (norm(quad)^2)/m];
        %prods = [prods quadratic_product(lambda,comp_z,1)/n];
    end
   

    mean_q2 = sum(prods);
    mean_q3 = sum(prods - prods_2 + length(mode)/m);


    var_q1 = var(prods_3*m);
    var_q2 = var(prods*m);
    var_q3 = var(-prods*m + prods_2*m - length(mode));
%    prods*m
%    prods*m - prods_2*m + length(mode)

%     col_new = generate_cov_row(x,par_1);
%     circ_new = [col_new; col_new(end-1:-1:2)];
%     [c_oo,c_uo,c_uu] = circ_partition(circ_new,length(x));
%     cond_var = (c_uu-c_uo*inv(c_oo)*c_uo');      
%     cstar = inv(toeplitz(circ_1));
%     cstar_lb = cstar(length(x)+1:end,length(x)+1:end);    
%     new = cstar_lb*(cond_var);
%     trace_term = trace(new);
%     q3_prods = -1*(- prods*m + prods_2*m - length(mode))
%     q2_prods = prods*m
%     trace_term

    %sigma_sq_q = 1;
    %mu_q  = mean(mode_z); 
    %comp_z'*inv(toeplitz(circ))*comp_z
%   likelihood_exact = -0.5*(trace_term+quad_m+sum(log(lambda))) 
%   likelihood_q1 = -0.5*(sum(prods_3)+sum(log(lambda)));
%   likelihood_q2 = -0.5*(mean_q2+quad_m+sum(log(lambda)));
    likelihood_q3 = -0.5*(mean_q3+quad_m+sum(log(lambda)));
       
%     exact_like = -0.5*(trace_term+quad_m+sum(log(lambda))) 
% 
%     new_like = -0.5*(mean_new+quad_m+sum(log(lambda)))
%        
%     old_like = -0.5*(mean_old+quad_m+sum(log(lambda)))
%        
%     old_old_like = -0.5*(sum(prods_3)+sum(log(lambda)))

    q = likelihood_q3;
end

end