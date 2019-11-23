function [q g] = mm_bound3(x,Y,p,n)
    %v = p(1:n-1);
    %w = p(n:2*n-2);
    %a = p(2*n-1:3*n-3);
    %fprintf(length(a))
    %h = p(3*n-2:3*n-1);
   
    
    v = p(1:n-1);
    w = p(n:2*n-2);
    a = p(2*n-1);
    h = p(2*n:2*n+1);

    
    if (h(1)<=0)||(h(2) <=0)
    q = -Inf;
    g = zeros(1,2*n);
    return
       
    else
    
    %eigenvalues
    
    col = generate_cov_row(x,h);
    circ = [col; col(end-1:-1:2)];
    lambda = ifft(circ);
    
    
    comp_w = [Y;w'];
    prod_w = quadratic_product(lambda,comp_w,1);
    %comp_w'*inv(toeplitz(circ))*comp_w
    
    comp_v = [zeros(n+1,1);v'];
    %comp_v
    %comp_v'*inv(toeplitz(circ))*comp_v
    prod_v = quadratic_product(lambda,comp_v,1);
    
    comp_trace = [zeros(2*n-1,1);1];
    %comp_trace'*inv(toeplitz(circ))*comp_trace
    prod_trace = quadratic_product(lambda,comp_trace,1);
    
    trace_term = prod_trace*(n-1); 
    
   
    q = prod_w+sum(log(lambda*2*n))+a*trace_term+prod_v-log(1+(norm(v)^2)/a)-(n-1)*log(a);
    if nargout > 1 % gradient required
        g = mm_bound3_g(x,Y,p,n);
        
    end
    
    %q = prod_w+sum(log(lambda))+trace_term+prod_v-log(1+sum(v.*(v./a)))-sum(log(a));
    %q=q-sum(log(h))*0.001;
   
end