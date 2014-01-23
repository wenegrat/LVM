function p = loadParams(ds)
% p.method = 5;

% Wind Data
% p.tlat = ds.ascat.lats;
% p.tlon = ds.ascat.lons(11:71); %40E-100E
% p.ttime = ds.ascat.time(:); 

p.tlat = ds.lat;
p.tlon = ds.lon;
p.ttime = ds.time;

% Model Parameters
p.lats = p.tlat;
p.lons = p.tlon;
% p.time = interp1(1:length(p.ttime), p.ttime, 1:.25:length(p.ttime));
% p.time = interp1(1:length(p.ttime), p.ttime, 1:.5:length(p.ttime));
p.time = p.ttime;

% p.ws = loadWS(p, ds.ascat.wstru);
p.ws = ds.taux;
p.ws(~isfinite(p.ws)) = 0;

p.maxVermodes = 2;
p.maxMermodes = 6;

p.damp = 1./(12*(30*24*60*60)); %12 month damping timescale for first mode


% Apply Reflectivity
% Note that applying a reflectivity directly during the integration is not
% quite equivalent to applying during the construction of SSH or U fields
% due to the finite differencing. This approach ensures
% that the reflected waves are truly damped by the proper percentage.
p.eref= .85;
p.wref = .85;


% Model Constants
p.beta = 2.28e-11;
p.deg = 110000; %convert from degrees to m
p.grav = 9.806;

p.timetosec = 60*60; %Daily hour wind stress

p.rho = 1025;
p.H = 4000;

% Constants for code progress display
p.tdisps = floor(length(p.time)/4);
p.wdisps = floor(length(p.lons).*length(p.time)/10);
end