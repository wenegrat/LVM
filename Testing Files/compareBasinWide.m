function compareBasinWide(ssh, ds, p)

tstrt = 180; % Discard first 6 months of transients

[nlon nlat ntime] = size(ssh(:,:,tstrt:end));
[~ , ~, ntimeobs] = size(ds.aviso.ssh(:,:,tstrt:end));

sshobs = gridObs(p.lons,p.lats,ds.aviso);

sshmvec = reshape(ssh(:,:,tstrt:end), nlon*nlat, ntime)*100;
sshovec = reshape(sshobs(:,:,tstrt:end), nlon*nlat, ntimeobs);

if (ntime > ntimeobs)
    lvec = 1:(ntimeobs);
else
    lvec = 1:(ntime);
end

b = NaN * zeros([nlon*nlat 2]);

[bf af] = butter(7, 2/(7), 'low');

for j=1:(nlon*nlat)
   sigmamod(j) = nanstd(sshmvec(j,:));
   sigmaobs(j) = nanstd(sshovec(j,:));
   
   mvshrt = squeeze(sshmvec(j,lvec));
   ovshrt = squeeze(sshovec(j,lvec));
   cr(j) = corr(filtfilt(bf, af, mvshrt'), filtfilt(bf, af, ovshrt'), 'rows', 'pairwise');
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
pcolor(p.lons, p.lats, sigmaobsm'); shading interp
colorbar
caxis([0 15]);
ylim(yl);

subplot(2,2,2)
pcolor(p.lons, p.lats, sigmamodm'); shading interp
colorbar
caxis([0 15]);
ylim(yl);

subplot(2,2,3)
pcolor(p.lons, p.lats, crm'); shading interp
caxis([-1 1]);
colorbar
ylim(yl);

subplot(2,2,4)
pcolor(p.lons, p.lats, bm'); shading interp
caxis([-2 2]);
ylim(yl);

colorbar
end