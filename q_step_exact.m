function q = q_step_exact(x,Y,mode,par,m,h_hat)

%generating circ row
col = generate_cov_row(x,par);
circ = [col; col(end-1:-1:2)];

global mu_q;
global sigma_sq_q;
global trace_term;



if (par(1)<=0)||(par(2) <=0)
    q = -Inf;
    
else
    %eigenvalues
    lambda = ifft(circ); 
    %iterate through samples
    %mean(mean(Z))
 
    
%   mu_q = [Y; mean(Z')'];
    mode_z = [Y;mode];
    

    
    mu =0;
    
    quad_m = sqrt(lambda).^-1.*ifft(mode_z-mu); 
    quad_m = norm(quad_m)^2;
  

    
    col_new = generate_cov_row(x,h_hat);
    circ_new = [col_new; col_new(end-1:-1:2)];
    [c_oo,c_uo,c_uu] = circ_partition(circ_new,length(x));
    cond_var = (c_uu-c_uo*inv(c_oo)*c_uo');      
    cstar = inv(toeplitz(circ));
    cstar_lb = cstar(length(x)+1:end,length(x)+1:end);    
    new = cstar_lb*(cond_var);
    trace_term = trace(new);
%     
    %sigma_sq_q = 1;
    %mu_q  = mean(mode_z); 
    %comp_z'*inv(toeplitz(circ))*comp_z
    
    
%     q_1 = -0.5*(sum(prods)+quad_m+sum(log(lambda)));
%     
%     q = q_1 + (0.5*(sum(prods_2)-length(mode)));

% 

       
%       exact_like = -0.5*(trace_term+quad_m+sum(log(lambda))) 
% 
%       new_like = -0.5*(mean_new+quad_m+sum(log(lambda)))
%       
%       old_like = -0.5*(mean_old+quad_m+sum(log(lambda)))
%       
%       old_old_like = -0.5*(sum(prods_3)+sum(log(lambda)))

      q = -0.5*(trace_term+quad_m+sum(log(lambda)));
end

end