function u = constructU(q, p, c, pmodes, option )
%===================================================
% Function constructU
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
[ndepths, ~] = size(pmodes);

[kelv ross force refl rere] = switchTypes(option);
            

%Make the q vector depending on options
qvec = makeqvec(q, p, force, refl, rere);

%Preallocate for speed
uvec = zeros([ndepths length(p.lons).*length(p.time) length(p.lats) ]);

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
                utemp = phik(y)./sqrt(2).*qt; 
                uvec(:,:,y) = uvec(:,:,y) + pmodes(:,n)*utemp;
                end
            elseif ross     %Check to see if the Rossby wave should be added
                h = sqrt( m*(m+1)/(2*(2*m + 1))).*(phip1(y)./sqrt(m+1) - phin1(y)./sqrt(m));
                utemp = qt.*h;
                uvec(:,:,y) = uvec(:,:,y)+ pmodes(:,n)*utemp;
            end
        end
     end
end
 %Reshape and permute SSH vector back into dims of: Lat x Lon x Time 
 u = reshape(uvec, ndepths,  length(p.lons), length(p.time), length(p.lats));
 u = permute(u, [1 2 4 3]);
end


