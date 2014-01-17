function sshs  = resampleSSH(ssh, p, resampleTime)

[nlons nlats nts] = size(ssh);

svec = reshape(ssh, nlons.*nlats, nts);

sveco = NaN * zeros([nlons*nlats length(resampleTime)]);
for i=1:(nlons*nlats)
%     sveco(i,:) = interp1(p.time, svec(i,:), resampleTime);
      sveco(i,:) = BlockMean(svec(i,:),1,7*2);
      sveco(i,:) = sveco(i,:) - nanmean(sveco(i,:));
end
sshs = reshape(sveco, nlons, nlats, length(resampleTime));
end