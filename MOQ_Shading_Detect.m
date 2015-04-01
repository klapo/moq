function DATstr = MOQ_Shade_Detect(DATstr)
% PURPOSE: 
% Finds shaded time steps using mean and standard deviation of transmissivity
% for angular bins of solar angles. Assumes that shaded times behave
% differently than non-shaded times (i.e., only applicable to fairly sunny
% locations)
%
% SYNTAX:
% 	SHADE_IND = MOQ_Shade_Detect(DATstr,SWFLAG,FIGFLAG,varargin)
%
% INPUTS:
%   DATstr  = 1x1 matlab structure: Contains the following fields
%           EL      = Nx1 vector - Elevation angle
%           AZ      = Nx1 vector - Azimuth angle
%           SWdwn   = Nx1 vector - Downwelling shortwave [Wm-2]
%   		SWFLAG  = Nx7 matrix of QC flags (from SW_Obs_QC_v2.m)
%
% OUTPUTS:
%	DATstr.SWFLAG	= Nx8 vector - QC flags for passing QC. Adds column for shading:
% 			Column 8 - Shading
%				1 = shaded
%				0 = not-shaded

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
if isfield(DATstr,'SWFLAG')
    SWFLAG = DATstr.SWFLAG;
else
    error('No azimuth angle was found (fieldname AZ)')
end

% Check for consistency
if length(SWdwn) ~= length(EL) || length(SWdwn) ~= length(AZ) || length(SWdwn) ~= length(SWFLAG)
    error('All data vectors must be the same length')
end

%%%%%%%%%%
%% CODE %%
%%%%%%%%%%

%% Search zenith and azimuth angular space for shading
% Only applicable if we have more than 3 years of data
obs_len = length(DATstr.(s{n}).SWdwn);
dt = Get_dt(DATstr.(s{n}).t);
dt_per_year = 365/dt;
num_years = obs_len/dt_per_year;
if num_years < 3
	disp('Function requires 3 years of observations. Exiting...')
	return
else

%% Bin transmissivity
% SWQCFLAG - general pass fail
[~,Tr_bar,Tr_std,bincount] = SW_Obs_QC_PlottingTool_SolGeo(DATstr.(s{n}),~DATstr.(s{n}).SWQCFLAG(:,7),0);

%% 2D space search
% Find points in time series within discrete solar geometry grid.
dEL = .5;                                           % Elevation bin width
ELplot = [dEL/2:dEL:90-dEL/2];                      % Middle of each EL bin
dAZ = .5;                                           % Azimuth bin width
AZplot = [dAZ/2:dAZ:360-dAZ/2];                     % Middle of each AZ bin
[ind_EL,ind_AZ] = find(Tr_bar > .4 | Tr_std > .15); % Index within the 2 dimensional space
ind = find(Tr_bar > .4 | Tr_std > .15);             % Linear index (reshaped to vector as MATLAB is want to do)
ind_exp = [];                                       % Pre-allocate

% Limit the 2D search space to only observed values
ELmax = round(max(DATstr.(s{n}).EL)./dEL).*dEL;
ELmin = 0;
AZmax = round(max(DATstr.(s{n}).AZ)./dAZ).*dAZ;
AZmin = round(min(DATstr.(s{n}).AZ)./dAZ).*dAZ;

% Discretized solar geometry space
for k = 1:length(ind_EL)
	temp = find(DATstr.(s{n}).AZ > AZplot(ind_AZ(k)) - dAZ/2 ...
		& DATstr.(s{n}).AZ <= AZplot(ind_AZ(k)) + dAZ/2 ...
		& DATstr.(s{n}).EL > ELplot(ind_EL(k)) - dEL/2 ... 
		& DATstr.(s{n}).EL <= ELplot(ind_EL(k))+ dEL/2);
	ind_exp = [ind_exp; temp];
end

% Times not classified as exposed (not a member of exposed index)
% 0 = not classified as shaded, 1 = classified as shaded
exc_ind = ~ismember(1:length(DATstr.(s{n}).AZ),ind_exp)' & DATstr.(s{n}).EL > 0;
DATstr.(s{n}).SWQCFLAG = [DATstr.(s{n}).SWQCFLAG,logical(exc_ind)];
