function gssh = gridObs(lons, lats, aviso)

[XI YI] = meshgrid(lats, lons);
[X Y] = meshgrid(aviso.lats, aviso.longs);

gssh = NaN * zeros([length(lons) length(lats) length(aviso.ssh)]);

for i=1:length(aviso.ssh)
    gssh(:,:,i) = interp2(X,Y, squeeze(aviso.ssh(:,:,i)), XI, YI);
end


end