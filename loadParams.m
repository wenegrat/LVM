function p = loadParams()

% Basic Run Parameters
% ===============================
p.maxVermodes = 2;          % How many vertical modes to use
p.maxMermodes = 6;          % How many meridional modes to use

% p.wp
% Should be true unless the F field already exists for the data and 
% mode choices (ie. re-running the same projection).
p.wp = true;

% Wind Stress Data
% ===============================
ds = ncdfread('/Users/JacobWenegrat/Documents/LWM/Wforc/ECMWF_TAUX_IO.nc');
p.tlat = ds.lat;
p.tlon = ds.lon;
p.ttime = ds.time;
p.timetosec = 60*60; % Update depending on unit of time (datestr -> timetosec = 24*60*60)
p.resampleTime = 7*2; % For resampling in time of output, number of consecutive timesteps to average.


% Model Parameters
% ===============================
% Model uses same horizontal grid as forcing
% Option to interpret the wind stress in time 
%
p.lats = p.tlat;
p.lons = p.tlon;
[p.ws p.time] = interpStress(ds.taux, p, 1); % Set interpolation here to meet CFL criteria


%
% Calculate Vertical Modes
% ===============================
% Can specify them from known data, or use calcVertModes to calculate
% based on stratification data.

% [wmd pmd c] = calcVertModes(ds, p);
p.pmd = [4.93 4.926 461.5127; 4.16 4.1592 353.9108]'; % ndepths x nmodes (nmodes >= p.maxVermodes)
p.c = [2.5335 1.6257];                                % Baroclinic phase speed.
p.vertdepths = [0 15 100];


% Output Options
% ========================================
% p.coption sets the type of fields to use for output:
%                       1 - All
%                       2 - Kelvin Forced 
%                       3 - Kelvin Reflected 
%                       4 - Kelvin Reflected x 2
%                       5 - Rossby Forced
%                       6 - Rossby Reflected
%                       7 - Rossby Reflected x 2
% ===================================================
p.coption = 1;


% Apply Reflectivity
%=========================================
% Reflectivity is applied to the Q fields during creation of SSH/U/etc.
% fields, and essentially mimics leaky boundaries.
%
% This is separate from the reflection coefficients that are a part of the 
% boundary conditions themselves which partion reflected energy into
% different modes.
p.eref= .85;
p.wref = .85;


% Timestepping Method
% ===================
% Ref: Durran (2013) - Numerical Methods for Fluid Dynamics
% 1 = Lax-Wendroff
% 2 = Upwind
p.method = 2;

% Forcing Integration Method
% ===================
% Select the type of numerical integration to use:
% 1 = Summation (fast, and consistent with Nagura 2010).
% 2 = Trapezoidal (slower, but more accurate)
p.intmethod = 1;


% Model Constants
% ========================================
p.beta = 2.28e-11;
p.deg = 110000; % convert from degrees to m
p.grav = 9.806;
p.damp = 1./(12*(30*24*60*60)); %12 month damping timescale for first mode
p.rho = 1025;
p.H = 4000;     % Mean ocean depth

% Constants for code progress display
p.tdisps = floor(length(p.time)/4);
p.wdisps = floor(length(p.lons).*length(p.time)/4);
end