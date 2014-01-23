function psi = hermitetest(m, y)

H(1,:) = ones(size(y)); %H0
H(2,:) = 2.*y;          %H1

for nn = 3:m+1          %H2-Hm
    H(nn,:) = 2.*(y'.*H(nn-1,:) - (nn-2).*H(nn-2,:));
end

tmp1 = 1;
for i=1:m
    tmp1 = tmp1*2;
end

psi = exp( (-1/2).*y.^2).*H(m+1,:)'./(sqrt(tmp1*factorial(m)*sqrt(pi)));
end