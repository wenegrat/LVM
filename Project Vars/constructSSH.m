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
            
qvecrere = reshape(squeeze(q(4,:,:,:,:)),  p.maxVermodes, p.maxMermodes + 1, length(p.lons).*length(p.time));


%Make the q vector depending on options
qvec = zeros([p.maxVermodes p.maxMermodes+1 length(p.lons).*length(p.time)]);

if (force)
    % Reshape the q vectors           
    qvecforc = reshape(squeeze(q(2,:,:,:,:)),  p.maxVermodes, p.maxMermodes + 1, length(p.lons).*length(p.time));
    qvec = qvec + qvecforc;
end
if (refl)
    qvecrefl = reshape(squeeze(q(3,:,:,:,:)),  p.maxVermodes, p.maxMermodes + 1, length(p.lons).*length(p.time));
    qvecrefl(:,1,:) = qvecrefl(:, 1,:).*p.wref; % Apply western BC reflectivity to Kelvin mode
    qvecrefl(:,2:end,:) = qvecrefl(:,2:end,:).*p.eref; % Apply Eastern BC reflectivity to all Rossby modes;
    qvec = qvec + qvecrefl;
end
if (rere)
    qvecrere = qvecrere.*p.wref.*p.eref; % Apply both BCs to the reflect + reflect;
    qvec = qvec + qvecrere;
end

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


% Helper function to isolate the switch statement
function [kelv ross force refl rere] = switchTypes(option)

    switch option
        case 1  %All forcing
                kelv = true;
                ross = true;

                force = true;
                refl = true;
                rere = true;
        case 2  %Kelvin Forced
                kelv = true;
                ross = false;

                force = true;
                refl = false;
                rere = false;
        case 3  %Kelvin Reflected
                kelv = true;
                ross = false;

                force = false;
                refl = true;
                rere = false;
        case 4  %Kelvin Reflected twice
                kelv = true;
                ross = false;

                force = false;
                refl = false;
                rere = true;
        case 5  %Rossby Forced
                kelv = false;
                ross = true;

                force = true;
                refl = false;
                rere = false;
        case 6  %Rossby Reflected
                kelv = false;
                ross = true;

                force = false;
                refl = true;
                rere = false;
        case 7  %Rossby Reflected twice
                kelv = false;
                ross = true;

                force = false;
                refl = false;
                rere = true;
    end
end

