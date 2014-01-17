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
pmd = [4.93 0; 4.16 0]';
c = [2.53 1.5];

lambda = sqrt(c./(p.beta));

% Preallocate for speed
if(wp)
F = NaN * zeros([p.maxVermodes p.maxMermodes+1 length(p.lons) length(p.time)]);
end
q = NaN * F;

for n=1:p.maxVermodes
    disp('')
    disp(['============= Baroclinic Mode Number: ', num2str(n),' ==========']);
    disp('')
    
    disp('============= Calculate Wind Stress Projection =========');
    if wp
    F(n,:,:,:) = windstressproject(p.ws, p, lambda(n), pmd(1,n));  % N x M x Lon x Time
    end

    disp('============= Integrate q function =====================');
    q(n,:,:,:) = integrateLVM_LW(squeeze(F(n,:,:,:)), p, c(n), pmd(1,n));
%         q(n,:,:,:) = integrateLVM(squeeze(F(n,:,:,:)), p, cn, dampn);


end

disp('============= Construct SSH Field  =====================');
ssh = constructSSH(q, p, c, pmd);
sshs = resampleSSH(ssh,p, motmodel.day.*24);

disp('============= Construct U Field    =====================');
    
disp(['Elapsed Time: ', num2str(toc./60), ' min']);


%% Plots
 cl = [-.2 .2];
%  
% pause
figure
       
       for i=1:length(motmodel.day)
       clf
       subplot(1,2,1)
       pcolor(p.lons, p.lats, squeeze(sshs(:,:,i))'.*.85);
       shading interp;
       caxis(cl);
       ylim([-10 10]);
       colorbar
       title(['t = ', num2str(i), '/', num2str(length(p.time))]);
       
       subplot(1,2,2)
%        pcolor(p.lons, p.lats, squeeze(ssh(:,:,i))');
%         pcolor(ds.aviso.longs(22:82), ds.aviso.lats, double(squeeze(ds.aviso.ssh(22:82,:,i)/100))'); 
        pcolor(motmodel.lon, motmodel.lat, double(squeeze(motmodel.SSH(:,:,i))'));
       shading interp;
       caxis(cl);
       ylim([-10 10]);
       colorbar
       title(['t = ', num2str(i), '/', num2str(length(p.time))]);
       drawnow 
       
       end
  
       %%
 cl = [-.2 .2];
 subplot(1,2,1)
 pcolor(p.lons, p.ttime, squeeze(sshs(:,26,:))'); 
 shading interp
 datetick('y','mm-yy');
 ylim([p.time(1) p.time(end)]);
 colorbar;
 caxis(cl);
 grid on
  subplot(1,2,2)
 pcolor(ds.aviso.longs(22:82), ds.aviso.time, double(squeeze(ds.aviso.ssh(22:82,27,2:end)/100))'); 
 shading interp
 datetick('y', 'mm-yy');
  ylim([p.time(1) p.time(end)]);
colorbar;
 caxis(cl*1);

% %        ylim([-300 300])
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