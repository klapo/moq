% This function uses air temperature data to estimate rainfall vs snowfall
% fractions.  Based on United States Army Corps of Engineers (1956)
% 
%RELEASE NOTES
% Written by Mark Raleigh (mraleig1@uw.edu)
% 
%SYNTAX
% [Rfall, Sfall] = RvsS(T,P,RvS)
% 
%INPUTS
% 1) T - Lx1 array of air temperature data (deg C)
% 2) P - Lx1 array of incremental precipitation data (mm/timestep)
% 3) RvS = control of how rain and snow are divided.
%         Enter 0 if you want a single temperature (PXTEMP) that divides rain and snow
%         Enter 1 if you want a linear transition between 2 temperatures (PXTEMP1 and PXTEMP2)
%         Enter 2 if you want no adjustments made to the precipitation
%         data.  All precipitation is considered snowfall.
% 
%OUTPUTS
% 1) Rfall = Lx1 array of rainfall (mm)
% 2) Sfall = Lx1 array of snowfall (mm)
% 
%NOTE: Default temperatures are set in the code to divide rain and snow.
%Check these and make modifications as appropriate.

function [Rfall, Sfall] = RvsS(T,P,RvS, v1, v2, varargin)

%% Default Parameters   

if nargin==3
    %Used if RvS = 0
    PXTEMP = 1;             %Temperature dividing rain from snow, deg C - if temp is less than or equal to PXTEMP, all precip is snow.  Otherwise it is rain.
    
    %Used if RvS = 1
    PXTEMP1 = -1;           %Lower Limit Temperature dividing tranistion from snow, deg C - if temp is less than or equal to PXTEMP1, all precip is snow.  Otherwise it is mixed linearly.
    PXTEMP2 = 3;            %Upper Limit Temperature dividing rain from transition, deg C - if temp is greater than or equal to PXTEMP2, all precip is rain.  Otherwise it is mixed linearly.
elseif nargin==4
    PXTEMP = v1;
    PXTEMP1 = [];
    PXTEMP2 = [];
elseif nargin==5
    PXTEMP = [];
    PXTEMP1 = v1;
    PXTEMP2 = v2;
else
    error('Invalid number of inputs')    
end

if (numel(PXTEMP1)~=1 && isempty(PXTEMP1)==0) || (numel(PXTEMP2) ~=1 && isempty(PXTEMP2)==0) || (numel(PXTEMP)~=1 && isempty(PXTEMP)==0)
    error('Invalid size of PXTEMP values')
end

if PXTEMP1 == PXTEMP2
    PXTEMP = PXTEMP1;
    RvS = 0;
    
else
    if PXTEMP1 > PXTEMP2
        disp('Warning: PXTEMP1 was greater than PXTEMP2.  Reversing them now')
        PXTEMP1a = PXTEMP1;
        PXTEMP1 = PXTEMP2;
        PXTEMP2 = PXTEMP1a;
        clear PXTEMP1a
    end
end

%% Checks

if numel(T) ~= numel(P)
    numel(T)
    numel(P)
    error('Error - temperature and precip arrays must be the same size')
end

if size(T,1) ~= 1
    if size(T,2) ~= 1
        error('Error - T must be an array')
    end
end

if size(P,1) ~= 1
    if size(P,2) ~= 1
        error('Error - P must be an array')
    end
end

% if nanmax(isnan(T))==1
%     error('NaN value detected in temperature.  These are not allowed')
% end
% 
% if nanmax(isnan(P))==1
%     error('NaN value detected in precipitation.  These are not allowed')
% end



%% Code

L = length(T);

transitionx = [PXTEMP1 PXTEMP2];
transitiony = [1 0];

Rfall = zeros(L,1);
Sfall = zeros(L,1);
 
 
for i = 1:L
    T_air = T(i);                 % air temperature at this time step (deg C)
    precip = P(i);                      % precipitation at this time step (mm)

% Divide Rain and Snow
if RvS ==0
    if T_air <= PXTEMP            % then the air temperature is cold enough for snow to occur
        fracsnow = 1.0;
    else                                % then the air temperature is warm enough for rain
        fracsnow = 0.0;
    end
    
elseif RvS ==1
    if T_air <= PXTEMP1
        fracsnow =1.0;
    elseif T_air >= PXTEMP2
        fracsnow =0.0;
    else
        fracsnow = interp1(transitionx, transitiony, T_air);
    end
elseif RvS ==2
    fracsnow = 1.0;
else
    error('Invalid rain vs snow option')
end

if fracsnow>1 || fracsnow<0
    error('Invalid snow fraction')
end


fracrain=1-fracsnow;

Rfall(i,1) = fracrain * precip;
Sfall(i,1) = fracsnow * precip;

end

