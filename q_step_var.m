function q = q_step_var(x,Y,Z,par,m,mu,n,mode,par_prev)


%generating circ row
col = generate_cov_row(x,par);
circ = [col; col(end-1:-1:2)];
prods = []; 
prods2 = [];

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
%     c_inv = inv(toeplitz(circ));
    for i = Z
        comp_z = [Y;i];
        %comp_z'*inv(toeplitz(circ))*comp_z
        quad = sqrt(lambda).^-1.*ifft(comp_z-mode_z);
        prods = [prods (norm(quad)^2)];
        %prods = [prods quadratic_product(lambda,comp_z,1)/n];
    end
    

%     for i = Z
%         comp_z = [Y;i];
%         %comp_z'*inv(toeplitz(circ))*comp_z
%         quad = sqrt(lambda).^-1.*ifft(comp_z); 
%         prods2 = [prods2 (norm(quad)^2)];
%         %prods = [prods quadratic_product(lambda,comp_z,1)/n];
%     end
%     
%     
%     
%     
%     
%     col_new = generate_cov_row(x,par_prev);
%     circ_new = [col_new; col_new(end-1:-1:2)];
%     [c_oo,c_uo,c_uu] = circ_partition(circ_new,length(x));
%     cond_var = (c_uu-c_uo*inv(c_oo)*c_uo');
% %     
% %     
%     M = zeros(length(comp_z));
%     M(length(Y)+1:end,length(Y)+1:end) = cond_var;
% 
%     c_inv = inv(toeplitz(circ));
%     c_inv = (c_inv+c_inv')/2;
%     CC = c_inv*M;
% 
%     variance = [var(prods2) var(prods)]
%     diff_act = var(prods2) - var(prods)
%     diff = 4*mode_z'*CC*c_inv*mode_z
%     trace_term = 2*trace(CC*CC)
%     
    
    
    
%     sigma_sq_q = 1;
%     mu_q  = mean(mode_z);
    
    quad_m = sqrt(lambda).^-1.*ifft(mode_z-mu); 
    quad_m = norm(quad_m)^2;

    q = -0.5*(sum(prods)/m+quad_m+sum(log(lambda*2*n)));
end

end