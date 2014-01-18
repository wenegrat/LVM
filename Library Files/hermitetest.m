function psi = hermitetest(m, y)

H(1,:) = ones(size(y));
H(2,:) = 2.*y;

for nn = 2:m
    H(nn+1,:) = 2.*(y.*H(nn+1-1) - (nn+1-1).*H(nn+1-2));
end

tmp1 = 1;
for i=1:m
    tmp1 = tmp1*2;
end

psi = exp( (-1/2).*y.^2)*H(m+1)/(sqrt(tmp1*factorial(m)*sqrt(pi)));
end