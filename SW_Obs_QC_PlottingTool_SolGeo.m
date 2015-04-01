function [h,Tr_bar,Tr_std,bincount] = SW_Obs_QC_PlottingTool_SolGeo(DATstr,SWFLAG,FIGFLAG,varargin)
% Graphs downwelling shortwave accoring to azimuth and solar zenith angle.
% Used to identify visually/subjectively potentially shaded solar geometery
% configurations.
%
% SYNTAX:
%	[h,Tr_bar,Tr_std,bincount] = SW_Obs_QC_PlottingTool_SolGeo(DATstr,SWFLAG)
%	[h,Tr_bar,Tr_std,bincount] = SW_Obs_QC_PlottingTool_SolGeo(DATstr,SWFLAG,FIGFLAG)
%
% INPUTS:
%	DATstr 	= 1x1 matlab structure: Contains the following fields
%			EL 		= Nx1 vector - Elevation angle
%			AZ 		= Nx1 vector - Azimuth angle
%			SWdwn 	= Nx1 vector - Downwelling shortwave [Wm-2]
%	SWFLAG 	= Nx7 matrix of QC flags (from SW_Obs_QC_v2.m)
%	FIGFLAG = Optional input, produce figures of binned mean and std tranmissivity?
%				1 = Yes, please, figured
%				0 = No (default)
%
% OUTPUTS:
%	h 		= 2x1 vector of figure handles. Returns NaNs if FIGFLAG = 0
%	Tr_bar 	= JxK matrix of mean transmissivity for each angular bin in the
%				space defined by the solar zenith and azimuth angles 
%	Tr_std 	= JxK matrix of transmissivity standard deviation for each angular
%				bin in the space defined by the solar zenith and azimuth angles
%	bincount= JxK matrix, number of time steps within given angular bin

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
if size(SWFLAG,2) > 1
	% If SWFLAG a matrix, collapse down to a good, no good index
	SWFLAG = sum(SWFLAG,2);
end

% Check for consistency
if length(SWdwn) ~= length(EL) || length(SWdwn) ~= length(AZ) || length(SWdwn) ~= length(SWFLAG)
	error('All data vectors must be the same length')
end

if nargin == 2
	FIGFLAG = 1;
end

%%%%%%%%%%%
%% GRAPH %%
%%%%%%%%%%%

% Transmissivity
dEL = .5;											% Elevation bin width
ELdiscrete = [0:dEL:90];							% Discrete elevation bin edges
ELplot = [dEL/2:dEL:90-dEL/2];						% Middle of each EL bin
dAZ = .5;											% Azimuth bin width
AZdiscrete = [0:dAZ:360];							% Discrete azimuth bin edges
AZplot = [dAZ/2:dAZ:360-dAZ/2];						% Middle of each AZ bin
TOA = 1365.*sind(EL);								% Extraterrestrial solar
Tr = SWdwn./TOA;									% Transmissivity
ind = find(EL > 0 & SWFLAG == 0 & ~isnan(SWdwn) & Tr > 0 & Tr < 1);	% Only plot day time points that pass QC
[Tr_bar,Tr_std,bincount] = bindata2(Tr(ind),EL(ind),AZ(ind),ELdiscrete,AZdiscrete); % Binned transmissiity according to gridded AZ and EL

% Remove Nan
[i,j] = find(~isnan(Tr_bar));
ind = find(~isnan(Tr_bar));

% Figures of binned tranmissivity
if FIGFLAG
	% Plot mean vs solar geo
	h(1) = figure;										% Handle to figure
	scatter(AZplot(j),90-ELplot(i),20,Tr_bar(ind),'filled')

	% Formating
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
	SetFigureProperties([1.5,1],gcf)

	% Plot standard deviation vs solar geo
	h(2) = figure;										% Handle to figure
	scatter(AZplot(j),90-ELplot(i),20,Tr_std(ind),'filled')

	% Formating
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
	SetFigureProperties([1.5,1],gcf)
else
	h = NaN(2,1);
end

