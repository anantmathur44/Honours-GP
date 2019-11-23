function q = mm_bound2(x,Y,p,n)
    v = p(1:n+1);
    w = p(n+2:2*n);
    a = p(2*n+1);
    h = p(2*n+1:2*n+2);
    
    if (h(1)<=0)||(h(2) <=0)
    q = -Inf;
    return
       
    else
    %eigenvalues
    
    col = generate_cov_row(x,h);
    circ = [col; col(end-1:-1:2)];
    lambda = ifft(circ);
    
    
    comp_w = [Y;w'];
    prod_w = quadratic_product(lambda,comp_w,1);
    %comp_w'*inv(toeplitz(circ))*comp_w
    
    comp_v = [v';zeros(n-1,1)];
    %comp_v'*inv(toeplitz(circ))*comp_v
    prod_v = quadratic_product(lambda,comp_v,0);
    

    trace_term = circ(1)*(n+1);
    
    
    
    


    q = prod_w+sum(log(lambda*2*n))+a*trace_term+prod_v-log(1+(norm(v)^2)/a)-(n+1)*log(a);
  
end