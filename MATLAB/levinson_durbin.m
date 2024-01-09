function f=levinson_durbin(gam,r)

% gam(1)*[1,rho(1:end)'], rho=gam(2:end)/gam(1);
% is the first row of the autocovariance matrix G,
% using Durbin's algorithm.

n=length(gam);
rho=gam(2:end)/gam(1);
r = r/gam(1);
phi = zeros(n,1);
psi = zeros(n,1);
phi(1) = r(1) ; psi(1) = rho(1);
for k=1:n-1
    phi_kr = phi(k:-1:1);
    psi_kr = psi(k:-1:1);
    rho_k = rho(1:k);  
    b = 1-rho_k'*psi(1:k);
    a_1 = (r(k+1)-rho_k'*phi_kr)/b;
    if k<n-1
        a_2 = (rho(k+1)-rho_k'*psi_kr)/b;
        psi(1:k)=psi(1:k)-a_2*psi_kr;
        psi(k+1)=a_2;
    end
    phi(1:k)=phi(1:k)-a_1*psi_kr;
    phi(k+1)=a_1;
end
f = phi;
