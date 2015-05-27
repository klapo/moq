function [SWFLAG,h] = SW_Obs_QC_v2(DATstr,FIG_FLAG,inverse_flag,varargin)
% Assigns QC flags for solar irradiance observations based on Long and
% Shi (2008) using the most stringent parameters given ("1st level" in the
% paper). While these may not be entirely appropriate for the
% climatological conditions of every location, calibrating the
% method to each site is generally not viable given the poor maintenance/
% conditions that are characteristic of mountain study areas. In addition to
% the Long and Shi method, snow burial is detected using  proxy meteorological
% variables and shading is detected using a transmissivity based detection.
%
% CREATED:
%	K. Lapo 06/12
% 
% MODIFICATION HISTORY:
%	Expanded algorithm, made more robust, specific to SW - K. Lapo 08/13	           
%	Added stand alone functiosn for detecting shading and snow burial
% 
% SYNTAX:
%	FLAGS = SW_Obs_QC_v2(DATstr)
%	FLAGS = SW_Obs_QC_v2(DATstr,FIG_FLAG)
%
% INPUTS:e
%	DATstr	= 1x1 structure 	- Contains fields with observations for QC
%		t		= Nx7 matrix, time_builder formatted dates
%		SWdwn	= Nx1 vector, downwelling irradiance (Wm^-2)
%		EL		= Nx1 vector, elevation angle (from the average cosine of
%					sza over the time step).
%		AZ		= Nx1 vector, solar azimuth angle
%		HA		= Nx1 vector, hour angle
%		SOLDIST = Nx1 vector, normalized Earth-Sun distance
%		T		= Nx1 vector, air temperature (C)
%		P		= Nx1 vector, precipitation during the time step (mm)
%		WIND	= Nx1 vector, average wind speed (m/s)
%		REF		= string,	string indicating how the data are referenced
%					to the time stamp: 'BEG','MID','END'
%		MF	 	= Kx1 vector, indices in time series of when maintenance was performed (optional)
%	FIG_FLAG= 1x1 logical, (default = 0) make QC figures%
%
% OUTPUTS:
%	FLAGS	= Nx8 vector - QC flags
%				Column 1 - Min global shortwave allowed
%				Column 2 - Night time check (shortwave)
%				Column 3 - Climatologically derived maximum shortwave (defunct, included only as a later reference) 
%				Column 4 - "Physical" maximum shortwave
%				Column 5 - Unphysically low transmissivity (shortwave)
%				Column 6 - Snow burial of instruments (solar and thermal)
%				Column 7 - Shading (requires 3 years of data) 
%				Column 8 - General pass (0) or fail (1). Logical reversed with 'inverse_flag' 
%	h		= 4x1 handle - Scatter plots of transmissivity in solar geometry space 

% inverse_flag makes the general pass/fail QC check have the opposite sign
% meaning as the other flags. 1 = passed QC, 0 = failed ANY test. Included for
% the sake of legacy uses of code

%%%%%%%%%%%%
%% CHECKS %%
%%%%%%%%%%%%
% Check if all required fields exist
if ~isfield(DATstr,'t') && size(DATstr.t,2) ~= 7
	error('Time matrix in time_builder format can not be found.')
else
	t = DATstr.t;
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
if isfield(DATstr,'T')
    T = DATstr.T;
else
	error('No air temperature was found (fieldname T)')
end
if isfield(DATstr,'P')
    P = DATstr.P;
else
	error('No precipitation was found (fieldname P)')
end
if isfield(DATstr,'WIND')
    WIND = DATstr.WIND;
else
	error('No wind was found (fieldname WIND)')
end
if isfield(DATstr,'REF')
	REF = DATstr.REF;
else
	error('No time stamp reference argument was found (fieldname REF)')
end
if isfield(DATstr,'MF')
	MF = DATstr.MF;
else
	MF = NaN;
end

% Check for consistency
if length(SWdwn) ~= length(EL) || length(SWdwn) ~= length(AZ) || length(SWdwn) ~= length(SOLDIST) || length(SWdwn) ~= length(HA)
	error('All data vectors must be the same length')
end
if length(SWdwn) ~= length(T) || length(SWdwn) ~= length(P) || length(SWdwn) ~= length(WIND)
	error('All data vectors must be the same length')
end

if nargin == 1
	FIG_FLAG = 0;									% Do not create QC figures
end
if nargin < 3
	inverse_flag = 0;
end

%%%%%%%%%%%%%%%%
%% Parameters %%
%%%%%%%%%%%%%%%%
% Physical parameters
sigma = 5.67*10^-8;									% Boltzmann's Constant
mew = sind(EL);										% cos(solar zenith angle) or sin(solar elevation angle)
S = 1367; 											% Solar Constant (at mean Earth-Sun distance) (W/m^2)
S = S .* SOLDIST.^2;								% Correct for solar distance
TOA = S.*mew;										% Solar flux @ TOA

% Adjustable parameters (most strict from Long and Shi 2008)
C1 = .92;
C4 = .95;
C5 = 100;
C6 = 465;
C7 = 120;
C8 = 590;
C11 = .65;
C13 = 14; % NOTE: Error in Long and Shi with this parameter?
C14 = 14; % NOTE: Error in Long and Shi with this parameter?
C15 = 200;
C16 = 25;

%%%%%%%%%%%%%%%
%% Algorithm %%
%%%%%%%%%%%%%%%
% pre-allocation
SWFLAG = zeros(length(SWdwn),7);

%% Climatological tests
% Long and Shi 2008 tests
IRind = SWdwn < -4;							% min global SW (if no IR leaking)
SWFLAG(IRind,1) = 1;
PNind = SWdwn > 10 & (mew == 0);			% check for night time
SWFLAG(PNind,2) = 1;
PMind = SWdwn > S .* 1.5 .* mew.^1.2 + 100;	% physical max global SW (hard maximum)
SWFLAG(PMind,4) = 1;					
% Defunct (Climate dependent max global SW) - regularly fails in winter at high elevaiton sites
CMind = SWdwn > S .* C1 .* mew.^1.2 + 50;	% climate dependent max global SW (soft maximum)
SWFLAG(CMind,3) = 1;
% Sub-Rayleigh scattering - THIS FIELD NEEDS TO BE CHANGED
SRind = SWdwn./TOA < 1/100 & EL > 15;
SWFLAG(SRind,5) = 1;						% Unphysically low values for higher solar elevation angles

%% Snow on dome
PARAMS.cumPcrit = 1;                            % Critical accumulated precip amount in interval (mm) 
PARAMS.int_len = 6/24;                          % Length of precip interval (serial days)
PARAMS.Tmelt = 1.5;                             % Snow melt threshold (degrees C)
PARAMS.Twind = -3;                              % Threshold between wet and sticky and dry and not (degrees C)
PARAMS.Wcrit_W = 3;                             % Wind threshold for wet snow (m/s)
PARAMS.Wcrit_C = 2;                             % Wind threshold for dry snow (m/s)
SWFLAG(:,6) = MOQ_SOD_Detect(DATstr,PARAMS);	% Snow-on-dome flag

%% Convert to logical
SWFLAG = logical(SWFLAG);						% Convert all flags to logical value (1 = did not pass given test, 0 = did pass)

%% Shading
% General QC pass for handing to shading function
SWFLAG(:,7) = logical(sum(SWFLAG,2));			% QC logical (0 = passed all, 1 = failed test(s))
[shadeflag,h] = MOQ_Shade_Detect(DATstr,FIG_FLAG);

%% General QC pass - Finalize for output 
SWFLAG = [SWFLAG(1:6),shadeflag];				% Cat all QC tests
if inverse_flag
	SWFLAG(:,8) = ~logical(sum(SWFLAG,2));		% General pass-fail (1 = passed all, 0 = failed test(s))
else
	% Default output
	SWFLAG(:,8) = logical(sum(SWFLAG,2));		% General pass-fail (1 = failed, 0 = passed)
end
