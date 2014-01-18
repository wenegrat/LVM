function psia = hermFuncAve(m, y)

dx1 = 2*pi/(2*m+1)/10;
nx1 = (y(end)-y(1))./dx1 + 1;

x1 = y(1)+(0:floor(nx1)-1)*dx1;

psi = hermiteeq(m, x1);

dx1exp = dx1.*ones(size(x1));

dxa = (y(3:end) - y(1:end-2))/2;
dxa = [dxa(1); dxa; dxa(end)];

psia = zeros(size(y));

for i=1:length(dxa)
    mask = (x1>=(y(i)-dxa(i)/2)) & (x1<(y(i)+dxa(i)/2));
    if sum(mask)~=0
    psia(i) = nansum(psi(mask).*dx1)./nansum(dx1exp(mask));
    end
end


end