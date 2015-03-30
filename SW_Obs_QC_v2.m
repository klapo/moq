function [SWFLAG,h1,h2] = SW_Obs_QC_v2(DATstr,SHADE_FLAG,TILT_FLAG)
% Assigns QC flags for solar irradiance observations based on Long and
% Shi (2008) using the most stringent parameters given ("1st level" in the
% paper). While these may not be entirely appropriate for the
% climatological conditions of every location, calibrating the
% method to each site is generally not viable given the poor maintenance/
% conditions that are characteristic of mountain study areas. In addition to
% the Long and Shi method, other control methods are employed based on 1) solar
% geometry and 2) proxy meteorological variables.
%
% CREATED:
%	K. Lapo 06/12
% 
% MODIFICATION HISTORY:
%	Expanded algorithm, made more robust, specific to SW - K. Lapo 08/13	           
% 
% SYNTAX:
%	FLAGS = SW_Obs_QC(DATstr)
%	FLAGS = SW_Obs_QC(DATstr,FIG_FLAG)
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
%		elev	= 1x1 scalar, site elevation asl (m)
%		REF		= string,	string indicating how the data are referenced
%					to the time stamp: 'BEG','MID','END'
%		MF	 	= Kx1 vector, indices in time series of when maintenance was performed
%	FIG_FLAG= 1x1 logical		- Indicates if QC figures are to be made
%
% OUTPUTS:
%	FLAGS	= Nx8 vector 		- QC flags
%	h1		= 1x1 handle		- Morn_v_Aft figure
%	h2		= 1x1 handle		- Hinkelman figure


%% Notes:
% This function needs to check for each flag provided and then the number
% of outputs being requested. If the flag is zero, but the handle is still
% being requested, output an empty value.

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
	% 	error('No maintenance index argument was found (fieldname MF)')
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
C13 = 14; % NOTE: Error in Long and Shi with this parameter.
C14 = 14; % NOTE: Error in Long and Shi with this parameter.
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
CMind = SWdwn > S .* C1 .* mew.^1.2 + 50;	% climate dependent max global SW (soft maximum)
SWFLAG(CMind,3) = 1;
PMind = SWdwn > S .* 1.5 .* mew.^1.2 + 100;	% physical max global SW (hard maximum)
SWFLAG(PMind,4) = 1;					
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

%% Tilting

%% General QC pass
SWFLAG(:,7) = ~logical(sum(SWFLAG,2));			% QC logical (1 = passed all, 0 = failed test(s))

%% Graphical tests for tilting and shading
% if FIG_FLAG
% 	h1 = SW_Obs_QC_PlottingTool_MornAft(DATstr,SWFLAG);% Morning vs afternoon irradiance
% 	h2 = SW_Obs_QC_PlottingTool_SolGeo(DATstr,SWFLAG);	% Solar geometry and irradiance
% else
% 	h1 = [];
% 	h2 = [];
% end
