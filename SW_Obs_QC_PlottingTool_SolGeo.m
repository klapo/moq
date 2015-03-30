function [h,Tr_bar,Tr_std,bincount] = SW_Obs_QC_PlottingTool_SolGeo(DATstr,SWFLAG)
% Graphs downwelling shortwave accoring to azimuth and solar zenith angle.
% Used to identify visually/subjectively potentially shaded solar geometery
% configurations.
%
% SYNTAX:
%	h = SW_Obs_QC_PlottingTool_SolGeo(DATstr,SWFLAG)
%
% INPUTS:
%
%
% OUTPUTS:
%

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
if size(SWFLAG,2) > 1
	% If SWFLAG a matrix, collapse down to a good, no good index
	SWFLAG = sum(SWFLAG,2);
end

% Check for consistency
if length(SWdwn) ~= length(EL) || length(SWdwn) ~= length(AZ) || length(SWdwn) ~= length(SOLDIST) || length(SWdwn) ~= length(HA) || length(SWdwn) ~= length(SWFLAG)
	error('All data vectors must be the same length')
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

% Old version of scatter plot -- very very slow and cumbersome
% scatter(AZ(ind),90-EL(ind),25,SWdwn(ind)./TOA,'filled')% Scatter plot of transmissivity vs solar geometry
% To plot in polar
% scatter((1-sind(EL(ind))).*sind(AZ(ind)),(1-sind(EL(ind))).*cosd(AZ(ind)),5,SWdwn(ind)./TOA,'filled')
% Old attempt at making a better figure for printing
% hm = mesh(AZplot,90-ELplot,Tr_bar,'mesh','column');
