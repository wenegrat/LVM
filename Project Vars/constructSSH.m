function ssh = constructSSH(q, p, c, pmodes, option )
%===================================================
% Function constructSSH
% Inputs:
%               q - from integrate function dim: 4 x N x M x Lat x Time
%               p - parameter struct
%               c - wave speeds
%               pmodes - psi array
%               Option for type of reconstruction:
%                       1 - All
%                       2 - Kelvin Forced 
%                       3 - Kelvin Reflected 
%                       4 - Kelvin Reflected x 2
%                       5 - Rossby Forced
%                       6 - Rossby Reflected
%                       7 - Rossby Reflected x 2
% ===================================================


[kelv ross force refl rere] = switchTypes(option);
            

%Make the q vector depending on options
qvec = makeqvec(q, p, force, refl, rere);

%Preallocate for speed
sshvec = zeros([length(p.lats) length(p.lons).*length(p.time)]);

for n=1:p.maxVermodes
     lambda = sqrt(c(n)./(p.beta));
     eta = p.lats*p.deg./lambda; % Non dimensionalized meridional dimension
     
    
     %Normalized Hermite Function
     for m=0:p.maxMermodes
         mm = m+1; % This is for indexing purposes
         if (m == 0)
            phik  = hermiteeq(0, eta); % Zeroeth hermite 
         else
             phip1 = hermiteeq(m+1, eta);
             phin1 = hermiteeq(m-1, eta);
         end
         qt = squeeze(qvec(n,mm,:))';
        for y = 1:length(eta)         
            if (m==0)
                if (kelv)   %Check to see if the Kelvin wave should be added
                sshvec(y,:) = sshvec(y,:) + (c(n)./p.grav)*phik(y)./sqrt(2).*qt.*pmodes(1,n); 
                end
            elseif ross     %Check to see if the Rossby wave should be added
                h = sqrt( m*(m+1)/(2*(2*m + 1))).*(phip1(y)./sqrt(m+1) + phin1(y)./sqrt(m));
                sshvec(y,:) = sshvec(y,:) + qt.*(c(n)./p.grav).*h.*pmodes(1,n);
            end
        end
     end
end
 %Reshape and permute SSH vector back into dims of: Lat x Lon x Time 
 ssh = reshape(sshvec, length(p.lats), length(p.lons), length(p.time));
 ssh = permute(ssh, [2 1 3]);
end


