function [sf,horiz] = SW_Obs_QC_Shading(DATstr,DEMpath,DEMname)
% Find angles (theta and phi) to horizon for shading applicaitons 
%
% SYNTAX:
%
%
% INPUTS:
%
%
% OUTPUTS:

%% Checks
if ~isdir(DEMpath)
	error('Path to DEM was not valid')
end

% Extract from data structure
if ~isfield(DATstr,'t') && size(DATstr.t,2) ~= 7
    error('Time matrix in time_builder format can not be found.')
else
    t = DATstr.t;
end
if isfield(DATstr,'lon')
    EL = DATstr.EL;
else
    error('No longitude was found (fieldname lon)')
end
if isfield(DATstr,'lat')
    EL = DATstr.EL;
else
    error('No longitude was found (fieldname lat)')
end
if isfield(DATstr,'SWdwn')
    SWdwn = DATstr.SWdwn;
else
    error('No downwelling irradiance was found (fieldname SWdwn)')
end
if isfield(DATstr,'EL')
    EL = DATstr.EL;
else
    error('No elevation angle was found (fieldname EL)')
end
if isfield(DATstr,'AZ')
    AZ = DATstr.AZ;
else
    error('No azimuth angle was found (fieldname AZ)')
end
if isfield(DATstr,'SOLDIST')
    SOLDIST = DATstr.SOLDIST;
else
    error('No normalized Earth-Sun distance was found (fieldname SOLDIST)')
end
if isfield(DATstr,'HA')
    HA = DATstr.HA;
else
    error('No solar hour angle was found (fieldname HA)')
end

% Check for consistency
if length(SWdwn) ~= length(EL) || length(SWdwn) ~= length(AZ) || length(SWdwn) ~= length(SOLDIST) || length(SWdwn) ~= length(HA)
    error('All data vectors must be the same length')
end
if length(SWdwn) ~= length(T) || length(SWdwn) ~= length(P) || length(SWdwn) ~= length(WIND)
    error('All data vectors must be the same length')
end

% Load DEM 
cd(DEMpath)
load(DEMname);
if 

	
% Find grid cell containing site in lat/long
[xm ym] = meshgrid(X,Y);							% Mesh grid (degrees)
d = ( (xm-lon).^2 + (ym-lat).^2).^(1./2);			% Distance from each grid center to site.
[r,c] = find(min(min(d)) == d);						% Row/Column of the grid containing site
[xm ym] = meshgrid(x,y);							% Mesh grid (UTM - meters)
xsite = xm(r,c);									% Easting of site grid box
ysite = ym(r,c);									% Northing of site grid box

% Calculate the slope to each grid cell
d = ( (xm-xm(r,c)).^2 + (ym-ym(r,c)).^2).^(1./2);	% Distance in meters
d(r,c) = 0;											% Force distance of site grid box to be zero
slope = (Z-Z(r,c))./d;								% Rise/run
theta = atand(slope);								% Convert slope to angle above the horizon

% For some discreet phi bin find the angle to the horizon
dphi = pi/36;										% 5 degree bins
ds = 25;											% Interpolated line with 25m spacing
distance = 0:ds:10000;								% Limit interpolated line to 10km
phi = pi/2:-dphi:-3*pi/2+dphi;						% Look from East clockwise
h = NaN(size(phi));
for k = 1:length(phi)
	xlin = distance.*cos(phi(k))+xsite;				% Project line from polar
	ylin = distance.*sin(phi(k))+ysite;				% "					"
	theta_proj = interp2(xm,ym,theta,xlin,ylin,'*nearest',NaN);	% Interpolate from polar line
	h(k) = nanmax(theta_proj);						% Maximum angle above the horizon is the visible horizon
end
h(h < 0) = 0;										% Remove angles below a flat horizon
