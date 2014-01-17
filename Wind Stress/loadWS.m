function out = loadWS(p, taux)

out = NaN * zeros([length(p.lons) length(p.tlat) length(p.time)]);

% Real Data
warning('off','all');
for i=1:length(p.tlat)
    for j=1:length(p.lons)
        out(j,i,:) = interp1(p.ttime, squeeze(taux(j+10,i,:)), p.time);
    end
end
warning('on','all');
% Fake Data
% for i=1:400
%     out(1:45,:,i) = .1*ones(size(out(1:45,:,i)));
% end

out(~isfinite(out)) = 0;
out  = double(out);
end