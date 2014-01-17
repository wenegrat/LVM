function compareBasinWideMODELED(ssh, mmod, p)

tstrt = 25; % Discard first 6 months of transients

[nlon nlat ntime] = size(ssh(:,:,tstrt:end));
% [~ , ~, ntimeobs] = size(ds.aviso.ssh(:,:,tstrt:end));


sshmvec = reshape(ssh(:,:,tstrt:end), nlon*nlat, ntime)*100;
sshovec = double(reshape(mmod(:,:,tstrt:end), nlon*nlat, ntime))*100;

lvec = 1:ntime;


b = NaN * zeros([nlon*nlat 2]);

[bf af] = butter(7, 2/(7), 'low');

for j=1:(nlon*nlat)
   sigmamod(j) = nanstd(sshmvec(j,:));
   sigmaobs(j) = nanstd(sshovec(j,:));
   
   mvshrt = squeeze(sshmvec(j,lvec));
   ovshrt = squeeze(sshovec(j,lvec));
%    cr(j) = corr(filtfilt(bf, af, mvshrt'), filtfilt(bf, af, ovshrt'), 'rows', 'pairwise');
      cr(j) = corr( mvshrt',  ovshrt', 'rows', 'pairwise');

   mask = isfinite(mvshrt) & isfinite(ovshrt);
   if (sum(mask) > 2)
   b(j,:) = regress(ovshrt(mask)', [mvshrt(mask)' ones(size(mvshrt(mask)'))]);
   end
end

sigmaobsm = reshape(sigmaobs, nlon, nlat);
sigmamodm = reshape(sigmamod, nlon, nlat);
crm = reshape(cr, nlon, nlat);
bm = reshape(b(:,1), nlon, nlat);


yl = [-15 15];

subplot(2,2,1)
contourf(p.lons, p.lats, sigmaobsm'); 
colorbar
caxis([0 15]);
ylim(yl);

subplot(2,2,2)
contourf(p.lons, p.lats, sigmamodm'); 
colorbar
caxis([0 15]);
ylim(yl);

subplot(2,2,3)
contourf(p.lons, p.lats, crm');
hold on
% contour(p.lons, p.lats, crm', [1 .8 .6], 'k');
hold off
caxis([0 1]);
colorbar
ylim(yl);

subplot(2,2,4)
contourf(p.lons, p.lats, bm', 0:.125:1.5);
caxis([0 1.5]);
ylim(yl);

colorbar
end