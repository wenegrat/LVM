function ssh = constructSSH(q, p, c, pmodes )


qvec = reshape(q,  p.maxVermodes, p.maxMermodes + 1, length(p.lons).*length(p.time));

%Preallocate for speed
sshvec = zeros([length(p.lats) length(p.lons).*length(p.time)]);

 for n=1:p.maxVermodes
     lambda = sqrt(c(n)./(p.beta));
     eta = p.lats*p.deg./lambda;
     
     phik  = hermiteeq(0, eta);
     %Normalized Hermite Function
     for m=0:p.maxMermodes
        mm = m+1;
        
     if (m > 0)
         phip1 = hermiteeq(m+1, eta);
         phin1 = hermiteeq(m-1, eta);
     end
        
        for y = 1:length(eta)
            if m ==0
                sshvec(y,:) = sshvec(y,:) + (c(n)./p.grav)*phik(y)./sqrt(2).*squeeze(qvec(n,mm,:))'.*pmodes(1,n); 
            else
                sshvec(y,:) = sshvec(y,:) + squeeze(qvec(n,mm,:))'.*(c(n)./p.grav)*( sqrt( m*(m+1)/(2*(2*m + 1))).*...
                    (phip1(y)./sqrt(m+1) + phin1(y)./sqrt(m+1))).*pmodes(1,n);

            end
        end
     end
     
     
 end
 ssh = reshape(sshvec, length(p.lats), length(p.lons), length(p.time));
 ssh = permute(ssh, [2 1 3]);
end