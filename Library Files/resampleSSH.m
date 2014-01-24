function [sshs resampletime] = resampleSSH(ssh, intime, resamplenum)

[nlons nlats nts] = size(ssh);

svec = reshape(ssh, nlons.*nlats, nts);

ntsre = floor(nts./resamplenum);

sveco = NaN * zeros([nlons*nlats ntsre]);
for i=1:(nlons*nlats)
      sveco(i,:) = BlockMean(svec(i,:),1,resamplenum); 
end
sshs = reshape(sveco, nlons, nlats, ntsre);

resampletime = intime(1:resamplenum:ntsre);

end