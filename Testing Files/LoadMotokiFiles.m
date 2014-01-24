motmodel = ncdfread('/Users/JacobWenegrat/Documents/LWM/Integration/output/SSH_out.nc');
%%
motq = ncdfread('/Users/JacobWenegrat/Documents/LWM/Integration/output/O_n1_f.nc');
%%
motwp = ncdfread('/Users/JacobWenegrat/Documents/LWM/Integration/output/I_n2.nc');
motin = ncdfread('/Users/JacobWenegrat/Documents/LWM/Wforc/ECMWF_TAUX_IO.nc');
%%
motu = ncdfread('/Users/JacobWenegrat/Documents/LWM/Integration/output/u15m_out.nc');
%%

for i=1:1682
   
    pcolor(motmodel.lon, motmodel.lat, double(squeeze(motmodel.SSH(:,:,i))'));
    shading interp;
    ylim([-15 15]);
    drawnow;
end
%%
% % mode = 1;
% qk  = squeeze(q(2,1,:,:,:));
% qk = permute(qk, [2 3 1]);
figure
xloc = 24;
subplot(1,3,1)
pcolor(squeeze(double(motmodel.SSH(xloc,:,:))));
shading interp
colorbar
caxis([-.15 .15]);

subplot(1,3,2)
pcolor(squeeze(sshs(xloc,:,:)));
shading interp;
colorbar
caxis([-.15 .15]);

yloc = 49;
subplot(1,3,3)
scatter(double(squeeze(motmodel.SSH(xloc,yloc,:))), squeeze(sshs(xloc,yloc,:)));
onetoone
title(num2str(regress(double(squeeze(motmodel.SSH(xloc,yloc,:))), squeeze(sshs(xloc, yloc,:)))))

%%
mode = 7;
qk  = squeeze(q(2,1,:,:,:));
qk = permute(qk, [2 3 1]);

subplot(1,3,1)
pcolor(squeeze(motq.coef(:,:,mode)));
shading interp
colorbar
caxis([-.15 .15]);

subplot(1,3,2)
pcolor(squeeze(qk(:,:,mode)));
shading interp;
colorbar
caxis([-.15 .15]);

xloc = 5;
tr = 4000:5000;
subplot(1,3,3)
scatter(double(squeeze(motq.coef(xloc,tr,mode))), squeeze(qk(xloc, tr,mode)));
onetoone
title(num2str(regress(double(squeeze(motq.coef(xloc,tr,mode)))', squeeze(qk(xloc, tr,mode))')))


%%

mode = 1;
Ik  = squeeze(F(2,:,:,:));
Ik = permute(Ik, [2 3 1]);

subplot(1,3,1)
pcolor(double(squeeze(motwp.wind_coef(:,:,mode)))*1e-3);
shading interp
colorbar
caxis([-.2 .2]);

subplot(1,3,2)
pcolor(squeeze(Ik(:,:,mode)));
shading interp;
colorbar
caxis([-.2 .2]);

xloc = 20;
subplot(1,3,3)
scatter(double(squeeze(motwp.wind_coef(xloc,1000:2000,mode))).*1e-3, squeeze(Ik(xloc, 1000:2000,mode)));
onetoone
title(num2str(regress(double(squeeze(motwp.wind_coef(xloc,1000:2000,mode)))'.*1e-3, squeeze(Ik(xloc, 1000:2000,mode))')))


%% Plots
 cl = [-.2 .2];
%  
% pause
figure
       
       for i=1:length(motmodel.day)
       clf
       subplot(1,2,1)
       pcolor(p.lons, p.lats, squeeze(sshs(:,:,i))');
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