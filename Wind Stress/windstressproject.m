function F = windstressproject(taux, p, lambda)

[nlons nlats ntimes] = size(taux);

yn = p.lats.*p.deg./lambda; % Normalized Meridional Dimension
dyn = yn(2) - yn(1);        % Delta y

tvec = double(reshape(permute(taux, [1 3 2]), nlons*ntimes, nlats));

%Preallocate
fvec = NaN * zeros([p.maxMermodes+1 nlons*ntimes]);
psivec = NaN * zeros([p.maxMermodes+1 length(yn)]);

% Make Herm Funcs
for m=0:p.maxMermodes+1
    psivec(m+1,:) = hermFuncAve(m, yn); %XXX possibly error in Mode 1...this follows Nagura
end
display('Hermite Polynomials Constructed');

%Loop over each time and x step.
for j=1:(nlons*ntimes)
    
    % Output Display
    if (mod(j, p.wdisps) == 0)
        display(['Stress Step: ', num2str(j), '/', num2str(nlons*ntimes)]);
    end
    
    for m=0:p.maxMermodes
       %Kelvin wave
       if m==0 

             psi = 1./sqrt(2).*psivec(m+1,:)';
             
            % Two options for type of numerical integration, see loadParams()
             if (p.intmethod == 1)
                fvec(m+1, j) = sum(psi.*squeeze(tvec(j,:))'.*dyn);
             elseif (p.intmethod == 2)
                fvec(m+1,j) =  trapz(yn, psi.*squeeze(tvec(j,:))');
             end
       %Rossby Wave     
       else    
           psin1 = psivec(m,:)';
           psip1 = psivec(m+2,:)';
           tint = sqrt( (m*(m+1))/(2*(2*m+1))).*...
               (psip1./(sqrt(m+1)) - psin1./sqrt(m)).*squeeze(tvec(j,:))';
           
             if (p.intmethod == 1)
                fvec(m+1, j) = sum(tint.*dyn);
             elseif (p.intmethod == 2)
                fvec(m+1,j) = trapz(yn, tint);
             end
       end
    end
end

F = reshape(fvec, p.maxMermodes+1, nlons, ntimes);

end