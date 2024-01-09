function [c_oo,c_uo,c_uu] = circ_partition(c,n)

C = toeplitz(c);
c_oo = C(1:n,1:n);
c_uo = C(n+1:end,1:n);
c_uu = C(n+1:end,n+1:end);
    
end