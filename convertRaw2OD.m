% dod = raw2OD( d )
%
% UI NAME:
% Intensity_to_OD 
%
% Converts internsity (raw data) to optical density
%
% INPUT
% d - intensity data (#time points x #data channels
%
% OUTPUT
% dod - the change in optical density

%Reference:
%to Homer2.
%

function dod = convertRaw2OD( d )

% convert to dod
dm = mean(abs(d),1);
nTpts = size(d,1);
dod = -log(abs(d)./(ones(nTpts,1)*dm));