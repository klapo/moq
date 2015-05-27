function [shadeflag,h] = MOQ_Shade_Detect(DATstr,FIG_FLAG,varargin)
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
%	FIG_FLAG = 1x1 logical, create figures (0 = default)
%
% OUTPUTS:
%	shadeflag	= Nx1 vector - QC Flag for shading to be added to QC matrix
%				1 = shaded
%				0 = not-shaded
%	h 		= 2x1 scalar, handles to plots of mean and standard deviation of
%			transmissivity. Returns NaN if FIG_FLAG is 0.

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

if nargin == 1
	FIG_FLAG = 0;
end

%%%%%%%%%%
%% CODE %%
%%%%%%%%%%

%% Search zenith and azimuth angular space for shading
% Only applicable if we have more than 3 years of data
obs_len = length(SWdwn);
dt = Get_dt(t);
dt_per_year = 365/dt;
num_years = obs_len/dt_per_year;
if num_years < 3
	disp('Shade detection requires 3 years of observations. Skipping...')
	h = NaN(4,1);
	shadeflag = zeros(length(SWFLAG),1);
	return
else

%% Bin transmissivity
% SWQCFLAG - general pass fail
if FIG_FLAG
	[h,Tr_bar,Tr_std,bincount] = SW_Obs_QC_PlottingTool_SolGeo(DATstr,DATstr.SWFLAG(:,7),FIG_FLAG);
else
	[~,Tr_bar,Tr_std,bincount] = SW_Obs_QC_PlottingTool_SolGeo(DATstr,DATstr.SWFLAG(:,7),0);
	h(1) = NaN; h(2) = NaN;
end

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
ELmax = round(max(EL)./dEL).*dEL;
ELmin = 0;
AZmax = round(max(AZ)./dAZ).*dAZ;
AZmin = round(min(AZ)./dAZ).*dAZ;

% Discretized solar geometry space
for k = 1:length(ind_EL)
	temp = find(AZ > AZplot(ind_AZ(k)) - dAZ/2 ...
		& AZ <= AZplot(ind_AZ(k)) + dAZ/2 ...
		& EL > ELplot(ind_EL(k)) - dEL/2 ... 
		& EL <= ELplot(ind_EL(k))+ dEL/2);
	ind_exp = [ind_exp; temp];
end

% Times not classified as exposed (not a member of exposed index)
% 0 = not classified as shaded, 1 = classified as shaded
exc_ind = ~ismember(1:length(AZ),ind_exp)' & EL > 0;
shadeflag = logical(exc_ind);

%% Build Figures 
if FIG_FLAG
    % Points that passed shading - mean transmissivity 
    h(3) = figure;
    scatter(AZplot(ind_AZ),90-ELplot(ind_EL),20,Tr_bar(ind),'filled')
    CLim = [0 1];
    set(gca,'CLim',CLim)
    COL = colorbar;
    cmap = cbrewer('seq','Reds',11);
    cmap = cmap(2:end,:);
    colormap(cmap)
    set(gca,'YDir','reverse')
    grid on
    xlabel('Azimuth')
    ylabel('SZA')
    ylabel(COL,'Transmissivity')
    axis([0 355 5 90])

    % Points that passed shading - standard deviation of transmissivity 
    h(4) = figure;
    scatter(AZplot(ind_AZ),90-ELplot(ind_EL),20,Tr_std(ind),'filled')
    CLim = [0 .3];
    set(gca,'CLim',CLim)
    COL = colorbar;
    cmap = cbrewer('seq','Reds',7);
    cmap = cmap(2:end,:);
    colormap(cmap)
    set(gca,'YDir','reverse')
    grid on
    xlabel('Azimuth')
    ylabel('SZA')
    ylabel(COL,'\sigma Transmissivity')
    axis([0 355 5 90])
else
    h(3) = NaN; h(4) = NaN; 
end
