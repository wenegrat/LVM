function us  = resampleU(u, resampleTime)

[ndepths nlons nlats nts] = size(u);

uvec = reshape(u, ndepths.*nlons.*nlats, nts);

uveco = NaN * zeros([ndepths.*nlons*nlats length(resampleTime)]);
for i=1:(ndepths.*nlons*nlats)

      uveco(i,:) = BlockMean(uvec(i,:),1,7*2); % For Model Comp

end
us = reshape(uveco, ndepths, nlons, nlats, length(resampleTime));
end