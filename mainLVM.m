%{
============== Linear Eq. Wave Model ====================
Solves the continuously stratified, wind-forced linear long wave model
on an Equatorial Beta Plane. (McCreary, 1981).

All user vars are in loadParams.

Output to outstruct.

Author: Jacob Wenegrat
wenegrat@uw.edu
Based on code from Motoki Nagura (2010).
%}
clc
close all
clear q

tic

p = loadParams();

lambda = sqrt(p.c./(p.beta));

% Preallocate for speed
if(p.wp)
F = NaN * zeros([p.maxVermodes p.maxMermodes+1 length(p.lons) length(p.time)]);
end
q = NaN * zeros([4 size(F)]);

disp(['Generate forcing field, set to: ', num2str(p.wp)])
disp('===========================================');
disp(['Number of Vertical Modes = ', num2str(p.maxVermodes)]);
disp(['Number of Meridional Modes = ', num2str(p.maxMermodes)]);
disp('===========================================');

for n=1:p.maxVermodes
    disp('')
    disp(['============= Baroclinic Mode Number: ', num2str(n),' ==========']);
    disp('')
    
    disp('============= Calculate Wind Stress Projection =========');
    if p.wp
    F(n,:,:,:) = windstressproject(p.ws, p, lambda(n));  % N x M x Lon x Time
    end

    disp('============= Integrate q function =====================');
    q(:,n,:,:,:) = integrateLVM_LW(squeeze(F(n,:,:,:)), p, p.c(n), p.pmd(1,n)); % 4 x N x M x Lon x Time
    % q(1...) is all (forced + all reflection)
    % q(2...) is forced only
    % q(3...) is single reflection
    % q(4...) refl+refl

    disp(['Time to calculate mode ', num2str(n),': ', num2str(toc./60), ' min']);
end

% Construct Dynamic Variables from q field.

disp('============= Construct SSH Field  =====================');
ssh = constructSSH(q, p, p.c, p.pmd, p.coption);
[sshs, ~] = resampleSSH(ssh, p.time, p.resampleTime);

disp(['Time to construct SSH ', num2str(n),': ', num2str(toc./60), ' min']);

disp('============= Construct U Field    =====================');
u = constructU(q, p, p.c, p.pmd, p.coption);
[us resampleTime] = resampleU(u, p.time, p.resampleTime);


% Create Output Struct
% ==========================
out.lats = p.lats;
out.lons = p.lons;
out.time = resampleTime;
out.depths = p.vertdepths;
out.Merid_modes = [-1 1:p.maxMermodes];
out.Vert_modes = 1:p.maxVermodes;
out.fproj = F;                      %Projection of forcing onto Merid. Modes
out.q = q;                          %Integrated q function
out.ssh = sshs;                     %ssh at original sampling frequency.
out.u   = us;                       %u at original sampling frequency.
out.params = p;
out.params = rmfield(out.params, 'ws'); %Remove forcing to keep size down.


disp(['Total Elapsed Time: ', num2str(toc./60), ' min']);



