function [tout timeout] = interpStress(taux, p, intnum)
% interpStress
% ============================
% Takes a forcing field and interprets it in time
%
% Inputs:
%           p = param struct
%           taux = Wind Stress in Lons x Lats x Time
%           intnum = Number of interpolating points ( > 1)
% Outputs:
%           tout = interpolated wind stress (in time)
%           timeout = interpolated time vector
% ============================

if (intnum == 1) % No interp
    tout = taux;
    timeout = p.ttime;
    return
end

% Interp
tout = NaN * zeros([length(p.lons) length(p.tlat) length(p.time).*intnum]);
timeout = interp1(1:length(p.ttime), p.ttime, 1:(1/intnum):length(p.ttime));

% XXX - Might be a performance gain to be had here. 
% warning('off','all');
for i=1:length(p.tlat)
    for j=1:length(p.lons)
        tout(j,i,:) = interp1(p.ttime, squeeze(taux(j,i,:)), timeout);
    end
end
% warning('on','all');

% Set areas over land etc. to zero wind stress.
tout(~isfinite(tout)) = 0;
tout  = double(tout);
end