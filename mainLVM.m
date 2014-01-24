%{
============== Linear Eq. Wave Model ====================
Based on the LVM Model of Motoki Nagura

%}
clc
close all
clear q

% Create Wind Stress Projection?
wp = false;
disp(['Generate forcing field, set to: ', num2str(wp)])

if (~exist('ds', 'var'))
        load('/Users/JacobWenegrat/Documents/IO/RawData/iodh.mat');
end

tic

p = loadParams(motin);

% Calculate Vertical Modes

% [wmd pmd c] = calcVertModes(ds, p);
pmd = [4.93 4.926; 4.16 4.1592]'; % in this form to match output of calcVertModes
c = [2.5335 1.6257];

lambda = sqrt(c./(p.beta));

% Preallocate for speed
if(wp)
F = NaN * zeros([p.maxVermodes p.maxMermodes+1 length(p.lons) length(p.time)]);
end
q = NaN * zeros([4 size(F)]);

disp('===========================================');
disp(['Number of Vertical Modes = ', num2str(p.maxVermodes)]);
disp(['Number of Meridional Modes = ', num2str(p.maxMermodes)]);
disp('===========================================');

for n=1:p.maxVermodes
    disp('')
    disp(['============= Baroclinic Mode Number: ', num2str(n),' ==========']);
    disp('')
    
    disp('============= Calculate Wind Stress Projection =========');
    if wp
    F(n,:,:,:) = windstressproject(p.ws, p, lambda(n));  % N x M x Lon x Time
    end

    disp('============= Integrate q function =====================');
    q(:,n,:,:,:) = integrateLVM_LW(squeeze(F(n,:,:,:)), p, c(n), pmd(1,n)); % 4 x N x M x Lon x Time
    % q(1... is all (forced and reflection)
    % q(2... is forced only
    % q(3... is single reflection
    % q(4...  refl+refl

    disp(['Time to calculate mode ', num2str(n),': ', num2str(toc./60), ' min']);
end

% Construct Dynamic Variables from q field.
%                       1 - All
%                       2 - Kelvin Forced 
%                       3 - Kelvin Reflected 
%                       4 - Kelvin Reflected x 2
%                       5 - Rossby Forced
%                       6 - Rossby Reflected
%                       7 - Rossby Reflected x 2
% ===================================================
coption = 1;

disp('============= Construct SSH Field  =====================');
ssh = constructSSH(q, p, c, pmd, coption); 
% sshs = resampleSSH(ssh,p, ds.aviso.time);
sshs = resampleSSH(ssh, p, motmodel.day*24);
    disp(['Time to construct SSH ', num2str(n),': ', num2str(toc./60), ' min']);

disp('============= Construct U Field    =====================');
u = constructU(q, p, c, pmd, coption);
us = resampleU(u, motmodel.day*24);

disp(['Total Elapsed Time: ', num2str(toc./60), ' min']);



