function ssh = constructSSH(q, p, c, pmodes )


qvecforc = reshape(squeeze(q(2,:,:,:,:)),  p.maxVermodes, p.maxMermodes + 1, length(p.lons).*length(p.time));
qvecrefl = reshape(squeeze(q(3,:,:,:,:)),  p.maxVermodes, p.maxMermodes + 1, length(p.lons).*length(p.time));
qvecrere = reshape(squeeze(q(4,:,:,:,:)),  p.maxVermodes, p.maxMermodes + 1, length(p.lons).*length(p.time));

qvecrefl(:,1,:) = qvecrefl(:, 1,:).*p.wref; % Apply western BC reflectivity;
qvecrefl(:,2:end,:) = qvecrefl(:,2:end,:).*p.eref; % Apply Eastern BC reflectivity to all Rossby modes;
qvecrere = qvecrere.*p.wref.*p.eref; % Apply both BCs to the reflect + reflect;

qvec = qvecforc;%XX + qvecrefl + qvecrere; % Total of forcing + reflected + 2x reflected.


%Preallocate for speed
sshvec = zeros([length(p.lats) length(p.lons).*length(p.time)]);


% psivec = NaN * zeros([p.maxMermodes+1 length(yn)]);
% for m=0:p.maxMermodes+1
%     psivec(m+1,:) = hermiteeq(m, yn);
% end

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
%          phip1 = phivec(m+2);
%          phin1 = phivec(m);
     end
        
        for y = 1:length(eta)
            
            % Note, here is where additional reconstructions (such as
            % reflected only) would take place.
            
            if m ==0
% XX                sshvec(y,:) = sshvec(y,:) + (c(n)./p.grav)*phik(y)./sqrt(2).*squeeze(qvec(n,mm,:))'.*pmodes(1,n); 
            else
                   h = sqrt( m*(m+1)/(2*(2*m + 1))).*(phip1(y)./sqrt(m+1) + phin1(y)./sqrt(m));
                   
                   sshvec(y,:) = sshvec(y,:) + squeeze(qvec(n,mm,:))'.*(c(n)./p.grav).*h.*pmodes(1,n);

            end
        end
     end
     
     
 end
 ssh = reshape(sshvec, length(p.lats), length(p.lons), length(p.time));
 ssh = permute(ssh, [2 1 3]);
end