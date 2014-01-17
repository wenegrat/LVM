function F = windstressproject(taux, p, lambda, psi0)
[nlons nlats ntimes] = size(taux);

% %Preallocate output
% F = NaN * zeros([p.maxMermodes+1 length(p.lons) length(p.time)]);
% 
% for t = 1:length(p.time)
%     
%     for x = 1:length(p.lons)
%         
%         for m= 1:p.maxMermodes+1
%             
%             if m==1
%                 % Kelvin Wave
%                 psi = hermiteeq(m-1, p.lats.*p.deg./lambda);
%                 F(m,x,t) = psi0./sqrt(2) * trapz(p.lats.*p.deg, psi.*squeeze(taux(x,:,t))'./lambda); % Kelvin Wave, YM99 Eq. 6a
%             else
%                 % Rossby Wave
%                 psin1 = hermiteeq(m-1, p.lats.*p.deg./lambda);
%                 psip1 = hermiteeq(m+1, p.lats.*p.deg./lambda);
%                 tint = (psip1./(sqrt(m+1)) - psin1./(sqrt(m))).*squeeze(taux(x,:,t))'./lambda;
%                 F(m,x,t) = psi0.*trapz(p.lats.*p.deg, tint);
%             end
%         end
%     end
% end


tvec = reshape(permute(taux, [1 3 2]), nlons*ntimes, nlats);
fvec = NaN * zeros([p.maxMermodes+1 nlons*ntimes]);
for j=1:(nlons*ntimes)
    for m=0:p.maxMermodes
       if m==0 %Kelvin wave
           psi = hermiteeq(m, p.lats.*p.deg./lambda);
           fvec(m+1,j) = 1./sqrt(2) * trapz(p.lats.*p.deg./lambda, psi.*squeeze(tvec(j,:))');
       else    %Rossby Wave
           psin1 = hermiteeq(m-1, p.lats.*p.deg./lambda);
           psip1 = hermiteeq(m+1, p.lats.*p.deg./lambda);
           tint = sqrt( (m*(m+1))/(2*2*m+1)).*...
               (psip1./(sqrt(m+1)) - psin1./sqrt(m)).*squeeze(tvec(j,:))';
           fvec(m+1,j) = trapz(p.lats.*p.deg./lambda, tint);
       end
    end
end

F = reshape(fvec, p.maxMermodes+1, nlons, ntimes);

end