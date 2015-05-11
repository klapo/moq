function [SF] = MOQ_SOD_Detect(DATstr,PARAMS)
% SOD module from SW_Obs_QC_v2.m
%
% SYNTAX:
%	SF = MOQ_SOD_Detect(DATstr)
%
% INPUTS:
%	DATstr	= 1x1 structure 	- Contains fields with observations for QC
%		t		= Nx7 matrix, time_builder formatted dates
%		SWdwn	= Nx1 vector, downwelling irradiance (Wm^-2)
%		T		= Nx1 vector, air temperature (C)
%		P		= Nx1 vector, precipitation during the time step (mm/hr)
%		WIND	= Nx1 vector, average wind speed (m/s)
%	PARAMS 	= 1x1 structure 	- Parameter values for SD algorithm
%
% OUTPUTS:
%	SF		= Nx1 vector, QC flags

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
% Check for consistency
if length(SWdwn) ~= length(T) || length(SWdwn) ~= length(P) || length(SWdwn) ~= length(WIND)
	error('All data vectors must be the same length')
end
% Is an index of maintenance being supplied?
if isfield(DATstr,'im')
	im = DATstr.im;
else
	im = [];
end

%%%%%%%%%%%%%%%%
%% Parameters %%
%%%%%%%%%%%%%%%%
dt_s = Get_dt(t);								% Time step in serial format
dlen = length(SWdwn);
SF = zeros(size(SWdwn));
[~,S] = RvsS(T,P,1);							% Determine amount of frozen precip

% Adjustable parameters 
cumPcrit = PARAMS.cumPcrit;						% Critical accumulated precip amount in interval (mm) 
int_len = PARAMS.int_len;						% Length of precip interval (serial days)
Tmelt = PARAMS.Tmelt;							% Snow melt threshold (degrees C)
Twind = PARAMS.Twind;							% Threshold between wet and sticky and dry and not (degrees C)
Wcrit_W = PARAMS.Wcrit_W;						% Wind threshold for wet snow (m/s)
Wcrit_C = PARAMS.Wcrit_C;						% Wind threshold for dry snow (m/s)

% Snow on dome
n_dt_day = round(int_len./dt_s);				% Number of time steps in a 6 hour period
PF = 0;											% Initialize precip flag
for n = 1+n_dt_day+1:dlen						% Loop over each observation point 
	ind = n-n_dt_day:n;							% Window of preceding 6 hours
	if intersect(ind,im)
		MF = intersect(ind,im);					% Maintenance performed during period in question
		ind = MF:ind(end);						
	else
		MF = 0;									% No maintenance performed
	end
	if ~PF										% Snow currently on the dome?
		PF = sum(S(ind)) > cumPcrit;			% Exceeded cumulative precip threshold over the interval 
	end
	if PF										% Precip event large enough to check
		TF = T(n) > Tmelt;		 				% Snow melt likely occured
		WF = (WIND(n) > Wcrit_W && T(n) >= Twind) ||...		% Wet snow blown off
			(WIND(n) > Wcrit_C && T(n) < Twind);			% Dry snow blown off
		if TF || WF || MF						% Snow melted, blown off, or maintenance performed
			PF = 0;								% No more snow on dome
		else
			SF(n) = 6;							% Flag for snow on dome 
		end
	end
end


