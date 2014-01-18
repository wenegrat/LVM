function h = hermiteeq(n,y)
%Note y should be non-dimensionalized by lambda

H0 = ones(size(y));
H1 = 2.*y;

psi0 = pi^(-.25).*exp(-(y.^2)/2).*H0;
psi1 = 2^(-1/2).*pi^(-.25).*exp(-(y.^2)/2).*H1;


if n == 0
    h= psi0;
elseif n == 1
    h = psi1;
else
    for m = 1:n-1
        psinp = y.*psi1 - sqrt(m/2).*psi0; % See Vallis AOFD 8.160b
        psinp = sqrt(2/(m+1)).*psinp;
        psi0 = psi1;
        psi1 = psinp;
    end
    h = psinp;
end