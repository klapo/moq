function hand = MOQ_LongShiLimits_Graph(DAT)
% Plots irradiances w/ QC limits, following Long and Shi 2008. 
% Reqires an input data structure and checks what irradiances are included.
%
% SYNTAX:
% 	hand = MOQ_LongShiLimits_Graph(DAT,outname)
%
% INPUTS:
%	DAT 	= 1x1 structure of met data, including irradiances. Requires data
%			structure to follow uniform naming structure used in other Lapo
%			code. 
%		Optional fields (must contain at least one of the following):
%			SWdwn 	= Nx1 vector of observed downwelling shortwave irradiance
%			SWup 	= Nx1 vector of observed reflected irradiances
%			LWdwn 	= Nx1 vector of observed downwelling longwave irradiance
%			LWup 	= Nx1 vector of observed upwelling longwave irradiance
%			
%		Required fields:
%			EL		= Nx1 vector of elevation angles (SW)
%			SOLDIST = Nx1 vector of normalized distance from the sun (SW)
%			T		= Nx1 vector of air temperature (preferentially from the
%					same height as the LW instrument)
%
% OUTPUTS:
%	hand 	= ?x1 vector, handles to plots generated in code

%% Parameters
% Long and Shi 2008: first level QC limits
C1 = .92; % NSA
C2 = .8; % NSA
C3 = .82; % SGP
C4 = .95; % SGP
C5 = 100;
C6 = 465;
C7 = 120;
C8 = 590;
C11 = .58;
C12 = 11;
C13 = 16;
C14 = 18;
C15 = 200;
C16 = 25;

% Physical parameters
sigma = 5.67*10^-8;

%% Flags/Checks
SWdwnFLAG = 0;
SWupFLAG = 0;
LWdwnFLAG = 0;
LWupFLAG = 0;

% SWdwn
if isfield(DAT,'SWdwn')
	SWdwnFLAG = 1;
	SWdwn = DAT.SWdwn;
	if isfield(DAT,'EL') && isfield(DAT,'SOLDIST')
		EL = DAT.EL;
		SZA = 90-EL;
		SOLDIST = DAT.SOLDIST;
	else
		error('Data structure does not contain EL and SOLDIST')
	end
end

% SWup
if isfield(DAT,'SWup') && SWdwnFLAG
	SWupFLAG = 1;
	SWup = DAT.SWup;
elseif isfield(DAT,'SWup') && ~SWdwnFLAG
	error('To plot SWup figure, SWdwn must be included')
end

% LWdwn
if isfield(DAT,'LWdwn') && sum(~isnan(DAT.LWdwn)) > 0
	LWdwnFLAG = 1;
	LWdwn = DAT.LWdwn;
	if isfield(DAT,'T')
		T = DAT.T+273.15; 
		Trange = linspace(min(T),max(T),100);
	else
		error('LWdwn must include T')
	end
end

% LWup
if isfield(DAT,'LWup') && sum(~isnan(DAT.LWup)) > 0
	LWupFLAG = 1;
	LWup = DAT.LWup;
	LWup_range = 150:450;
	if isfield(DAT,'T')
		T = DAT.T+273.15;
		Trange = 240:300;
	else
		error('LWup must include T')
	end
end

% Final check
if SWdwnFLAG + SWupFLAG + LWdwnFLAG + LWupFLAG == 0
	error('No irradiances found in data structure')
end

%% Figures
hand = NaN(5,1);

% SWdwn limits
if SWdwnFLAG 
    % Level 1 limits
    lvl1_lim =  1368 .* C1 .* sind(EL).^1.2 .* SOLDIST + 50;
    % Physically Possible limits
    phys_lim = 1368 .* sind(EL).^1.2 .* SOLDIST .* 1.5 +  100;
    
    % Basic Plotting
    hand(1) = figure;
    plot(SZA,SWdwn,'.')
    hold all
    plot(SZA,lvl1_lim)
    plot(SZA,phys_lim)
    legend('Observed SW \downarrow','Climatological limits','Physical limits')
    xlabel('Solar Zenith Angle')
    ylabel('SW \downarrow (Wm^{-2})')
    grid on
	axis([0, max(SZA), 0, 1400])

    % Find number of points exceeding limits
    ind_lvl1 = find(SWdwn > lvl1_lim);
    ind_phys = find(SWdwn > phys_lim);
    str1(1) = {['Physical limit: ',num2str(length(ind_phys)./length(SWdwn)*100,'%0.3g'),'%']};
    str1(2) = {['Climatological limit: ',num2str(length(ind_lvl1)./length(SWdwn)*100,'%0.3g'),'%']};
    text(50,1000,str1)

end 

% Ratio of SWup to SWdwn
if SWupFLAG && SWdwnFLAG
    hand(2) = figure;
    plot(SWdwn,SWup,'.')
    hold all
    irrad_lims = [0 1400];
    plot(irrad_lims,irrad_lims.*.9)
    plot(irrad_lims,irrad_lims.*.2)
    legend('obs','Snow Albedo (.9)','Ground Albedo (.2)')
    xlabel('SW \downarrow (Wm^{-2})')
    ylabel('SW \uparrow (Wm^{-2})')
	axis([0 1200 0 1000])
	grid on
end

% LWdwn and Tair
if LWdwnFLAG 
    % Limits 
    lvl1_lim_low = C11 .* sigma .* Trange.^4;
    lvl1_lim_up = C12 + sigma .* Trange.^4;
   
   	hand(3) = figure;
    plot(T,LWdwn,'.')
    hold all
    plot(Trange,lvl1_lim_low)
    plot(Trange,lvl1_lim_up)
    legend('LW \downarrow obs','Low Clim Lim','High Clim Lim','Location','NorthWest')
    xlabel('T_{air} (K)')
    ylabel('LW \downarrow (Wm^{-2})')
	axis([Trange(1) Trange(end) 50 450])
	grid on
	
	% Find number of points exceeding limits
	% Limits w/ obs T
    lvl1_lim_low = C11 .* sigma .* T.^4;
    lvl1_lim_up = C12 + sigma .* T.^4;
    ind_low = find(LWdwn < lvl1_lim_low);
    ind_up = find(LWdwn > lvl1_lim_up);
    str1(1) = {['Lower limit: ',num2str(length(ind_low)./length(LWdwn)*100,'%0.3g'),'%']};
    str1(2) = {['Upper limit: ',num2str(length(ind_up)./length(LWdwn)*100,'%0.3g'),'%']};
    text(280,100,str1) 
end

% LWup and Tair
if LWupFLAG 
    % Limits
    lvl1_lim_low = sigma.*(Trange - C13).^4;
    lvl1_lim_up = sigma.*(Trange + C14).^4;
    
    hand(4) = figure;
	plot(T,LWup,'.')
    hold all
    plot(Trange,lvl1_lim_low)
    plot(Trange,lvl1_lim_up)
    legend('LW \uparrow obs','Low Clim Lim','High Clim Lim','Location','NorthWest')
    xlabel('T_{air} (K)')
    ylabel('LW \uparrow')
	axis([Trange(1) Trange(end) 50 450])
	grid on
	
	% Find number of points exceeding limits
	% Limits w/ obs T
	lvl1_lim_low = sigma.*(T - C13).^4;
    lvl1_lim_up = sigma.*(T + C14).^4;
	ind_low = find(LWup < lvl1_lim_low);
    ind_up = find(LWup > lvl1_lim_up);
    str1(1) = {['Lower limit: ',num2str(length(ind_low)./length(LWup)*100,'%0.3g'),'%']};
    str1(2) = {['Upper limit: ',num2str(length(ind_up)./length(LWup)*100,'%0.3g'),'%']};
    text(280,100,str1) 
end

% LWdwn vs LWup
if LWupFLAG && LWdwnFLAG    
	% Limits
	lvl1_lim_low = LWup_range - C15;
	lvl1_lim_up = LWup_range + C16;

    hand(5) = figure;
    plot(LWup,LWdwn,'.')
	hold all
	plot(LWup_range,lvl1_lim_low)
	plot(LWup_range,lvl1_lim_up)
    legend('Obs','Low Clim Lim','High Clim Lim','Location','NorthWest')
    xlabel('LW \uparrow (Wm^{-2})')
    ylabel('LW \downarrow (Wm^{-2})')
	axis([LWup_range(1) LWup_range(end) 50 540])
	grid on
	
	% Find number of points exceeding limits
	% Limits w/ obs LWup
	lvl1_lim_low = LWup - C15;
	lvl1_lim_up = LWup + C16;
	ind_low = find(LWdwn < lvl1_lim_low);
    ind_up = find(LWdwn > lvl1_lim_up);
    str1(1) = {['Lower limit: ',num2str(length(ind_low)./length(LWdwn)*100,'%0.3g'),'%']};
    str1(2) = {['Upper limit: ',num2str(length(ind_up)./length(LWdwn)*100,'%0.3g'),'%']};
    text(350,100,str1) 
end 
