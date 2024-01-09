function g = mm_bound3_g(x,Y,p,n)
    %v = p(1:n-1);
    %w = p(n:2*n-2);
    %a = p(2*n-1:3*n-3);
    %fprintf(length(a))
    %h = p(3*n-2:3*n-1);
   
    
    v = p(1:n-1);
    w = p(n:2*n-2);
    a = p(2*n-1);
    h = p(2*n:2*n+1);
    
    col = generate_cov_row(x,h);
    circ = [col; col(end-1:-1:2)];
    lambda = ifft(circ);
    
    %GRADIENT v
    
    v_pad = [zeros(n+1,1);v'];
    v_mm = matrix_vector_multi(v_pad,(lambda * 2 * n).^-1);
    v_mm = v_mm(n+2:end);
    g(1:n-1) = 2*v_mm -((1+(norm(v)^2)/a)^(-1))*(2/a)*v';

    %GRADIENT w
    comp_w = [Y;w'];
    w_mm = matrix_vector_multi(comp_w,(lambda * 2 * n).^-1);
    w_mm = w_mm(n+2:end);
    
    g(n:2*n-2) = 2*w_mm;
    
    %GRADIENT a

    comp_trace = [zeros(2*n-1,1);1];
    prod_trace = quadratic_product(lambda,comp_trace,1);
    trace_term = prod_trace*(n-1);
    da = trace_term + (norm(v)^2)/(a^2+a*norm(v)^2)-(n-1)/a;
    g(2*n-1) = da;
    
    %GRADIENT h
    delta = sqrt(eps);
    f_1 = mm_bound4(x,Y,[v w a h(1) h(2)],n);
    f_2 = mm_bound4(x,Y,[v w a (h(1)+delta) h(2)],n);
    dh1 = (f_2-f_1)/delta;
    f_2 = mm_bound4(x,Y,[v w a h(1) (h(2)+delta)],n);
    dh2 = (f_2-f_1)/delta;
    g(2*n) = dh1;
    
   
end