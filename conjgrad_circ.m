function x = conjgrad_circ(lambda,b,n,max_iterations)


%make sure x is column

b_pad = [b; zeros(n-1,1)];
x = matrix_vector_multi(b_pad,(lambda * 2 * n).^-1);
x = x(1:n+1);

x_pad = [x ; zeros(n-1,1)];
Ax = matrix_vector_multi(x_pad,lambda * 2 * n);
r = b - Ax(1:n+1); 


r_pad = [r ; zeros(n-1,1)];
p = matrix_vector_multi(r_pad,(lambda * 2 * n).^-1); % matrix multiplication
p = p(1:n+1);




for i = 1:max_iterations
 
    r_pad = [r ; zeros(n-1,1)];
    p_pad = [p ; zeros(n-1,1)];
    
    quad_rMr = quadratic_product(lambda,r_pad,1);
    quad_pAp = quadratic_product(lambda,p_pad,0);
    
    alpha = quad_rMr / quad_pAp;
    x = x + alpha * p;
    
    Ap = matrix_vector_multi(p_pad,(lambda * 2 * n)); 
    
    r = r - alpha * Ap(1:n+1);
    r_pad = [r ; zeros(n-1,1)];
    
    rnorm = r' * r;  
    bnorm = b'* b;
    
    if sqrt(rnorm)/sqrt(bnorm) < 1e-10
       
        break;
    end
    
    b = quadratic_product(lambda,r_pad,1)/quad_rMr;
    Mr = matrix_vector_multi(r_pad,(lambda * 2 * n).^-1);
    p = Mr(1:n+1)  + b*p;
end



if sqrt(rnorm) > 1e-10
    fprintf("Max Iterations PCG");
end

