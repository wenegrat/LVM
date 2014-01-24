function [us resampletime] = resampleU(u, intime, resamplenum)

[ndepths nlons nlats nts] = size(u);

uvec = reshape(u, ndepths.*nlons.*nlats, nts);

ntsre = floor(nts./resamplenum);

uveco = NaN * zeros([ndepths.*nlons*nlats ntsre]);
for i=1:(ndepths.*nlons*nlats)

      uveco(i,:) = BlockMean(uvec(i,:),1,resamplenum); % For Model Comp

end
us = reshape(uveco, ndepths, nlons, nlats, ntsre);

resampletime = intime(1:resamplenum:ntsre);
end