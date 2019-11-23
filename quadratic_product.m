function q = quadratic_product(lambda,x,inverse)
        %comp_z'*inv(toeplitz(circ))*comp_z
        if inverse
            quad = sqrt(lambda).^-1.*ifft(x); 
        else 
            quad = sqrt(lambda).*fft(x); 
        end
        q =  (norm(quad)^2);
        
        
        
end