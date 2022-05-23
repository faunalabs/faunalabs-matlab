%% epsilon · molar extinction coefficient conversion from known ua

%   Dave Haas, 6 July 2021

clc;

%% try this formula

%   eta = ( mua [1/cm] * 64500 [g/M] * 0.001 [M / mM] ) / (2.303 * Hb [g / L]

%   use some assumptions here, e.g.:

Hb_gPerLiter = 140;     % 

Hb_gPerM = 64450;       % Cope 1991 disseration appendix B
%Hb_gPerM = 64500;      % Prahl 1998
%Hb_gPerM = 66800;       % Kim et al 2005

MpermM = 0.001;

mua_HbO2_660 = 0.15;    % Bosschaart et al 2014 in units [mm ^ -1]
mua_HbR_660 = 1.64;     % Bosschaart et al 2014 in units [mm ^ -1]
mua_HbO2_995 = 0.65;    % Bosschaart et al 2014 in units [mm ^ -1]
mua_HbR_995 = 0.23;     % Bosschaart et al 2014 in units [mm ^ -1]

mua_mm2cm = 10;         % adjustment factor for mm^-1 to cm^-1 conversion

%muaConvFactor = 2.303;  % Prahl 1998, ua = 2.303 · epsilon · [C] · L
muaConvFactor = 2.3;    % Kim and Liu 2007, ua = 2.3 · epsilon · [C]
 
% these from Cope's Appendix B dissertation 1991
epsilonHbO2_660_expected = 0.3346;     % [mM ^ -1 · cm ^ -1]
epsilonHbR_660_expected = 3.4408;      % [mM ^ -1 · cm ^ -1]
epsilonHbO2_995_expected = 1.1848;   % [mM ^ -1 · cm ^ -1]
epsilonHbR_995_expected = 0.3047;    % [mM ^ -1 · cm ^ -1]

% lambda = 660 nm

epsilonHbO2_660_calc = ( (mua_HbO2_660 * mua_mm2cm) * Hb_gPerM * 0.001 ) / ...
    ( muaConvFactor * Hb_gPerLiter ) ;

epsilonHbR_660_calc = ( (mua_HbR_660 * mua_mm2cm) * Hb_gPerM * 0.001 ) / ...
    ( muaConvFactor * Hb_gPerLiter ) ;

% lambda = 995 nm

epsilonHbO2_995_calc = ( (mua_HbO2_995 * mua_mm2cm) * Hb_gPerM * 0.001 ) / ...
    ( 2.303 * Hb_gPerLiter ) ;

epsilonHbR_995_calc = ( (mua_HbR_995 * mua_mm2cm) * Hb_gPerM * 0.001 ) / ...
    ( 2.303 * Hb_gPerLiter ) ;

