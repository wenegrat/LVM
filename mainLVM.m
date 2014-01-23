%{
============== Linear Eq. Wave Model ====================
Based on the LVM Model of Motoki Nagura

%}
clc
close all
clear q

% Create Wind Stress Projection?
wp = true;
disp(['Generate forcing field, set to: ', num2str(wp)])

if (~exist('ds', 'var'))
        load('/Users/JacobWenegrat/Documents/IO/RawData/iodh.mat');
end

tic

p = loadParams(motin);

% Calculate Vertical Modes

% [wmd pmd c] = calcVertModes(ds, p);
pmd = [4.93 0; 4.16 0]'; % in this form to match output of calcVertModes
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

disp('============= Construct SSH Field  =====================');
ssh = constructSSH(q, p, c, pmd, 1); % Default is to construct full ssh field
% sshs = resampleSSH(ssh,p, ds.aviso.time);
sshs = resampleSSH(ssh, p, motmodel.day*24);

disp('============= Construct U Field    =====================');
    
disp(['Total Elapsed Time: ', num2str(toc./60), ' min']);


%% Plots
%  cl = [-.2 .2];
% %  
% % pause
% figure
%        
%        for i=1:length(motmodel.day)
%        clf
%        subplot(1,2,1)
%        pcolor(p.lons, p.lats, squeeze(sshs(:,:,i))');
%        shading interp;
%        caxis(cl);
%        ylim([-10 10]);
%        colorbar
%        title(['t = ', num2str(i), '/', num2str(length(p.time))]);
%        
%        subplot(1,2,2)
% %        pcolor(p.lons, p.lats, squeeze(ssh(:,:,i))');
% %         pcolor(ds.aviso.longs(22:82), ds.aviso.lats, double(squeeze(ds.aviso.ssh(22:82,:,i)/100))'); 
%         pcolor(motmodel.lon, motmodel.lat, double(squeeze(motmodel.SSH(:,:,i))'));
%        shading interp;
%        caxis(cl);
%        ylim([-10 10]);
%        colorbar
%        title(['t = ', num2str(i), '/', num2str(length(p.time))]);
%        drawnow 
%        
%        end
  
       %%
%  cl = [-.2 .2];
%  subplot(1,2,1)
%  pcolor(p.lons, p.ttime(1:513), squeeze(sshs(:,26,:))'); 
%  shading interp
%  datetick('y','mm-yy');
%  ylim([p.time(1) p.time(end)]);
%  colorbar;
%  caxis(cl);
%  grid on
%   subplot(1,2,2)
%  pcolor(ds.aviso.longs(22:82), ds.aviso.time, double(squeeze(ds.aviso.ssh(22:82,27,2:end)/100))'); 
%  shading interp
%  datetick('y', 'mm-yy');
%   ylim([p.time(1) p.time(end)]);
% colorbar;
%  caxis(cl*1);
% 
% % %        ylim([-300 300])
% %        xlim([x(1)./p.deg x(end)./p.deg]);
%        subplot(3,1,2:3)
%        plot(x./p.deg,u0,'b:',x./p.deg,ur(1,:),'r.-')
%        hold on
%               plot(x./p.deg,ur(2,:),'b.-')
%               plot(x./p.deg, ur(4,:), 'g');
%         hold off
%     %    axis([-1 1 -.1 1.1])
%        xlim([x(1)./p.deg x(end)./p.deg]);
%        ylim([-.5*1e8 .5*1e8]);
%        title(sprintf('%s , t=%0.2f','Lax-Wendroff', tn))
%        drawnow