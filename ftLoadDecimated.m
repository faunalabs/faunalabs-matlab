function [O, K, P, CUE, TAG, TRIAL] = ftLoadDecimated(tag)
%
%    [O, K, P, CUE, TAG, TRIAL] = ftLoadDecimated(tag);
%
%    Load the post-processed optical, kinematic, depth, and temperature
%    sensor data file associated with a tag deployment.
%
%    Dave Haas
%    dave.haas@faunalabs.com
%    21 July 2021

if nargin<1
   help ftLoadDecimated
   return
end

O=[]; K=[]; P=[]; CUE=[]; TAG=[]; TRIAL=[];

global TAG_PATHS;

if isempty(tag)
    fprintf('Specify a valid tag trial, e.g.: tt21_142a and try again.\n');
    return;
end

prhFileTarget = sprintf('%s/%sprh.mat', TAG_PATHS.PRH, tag);

% check if the file exists
if ~exist(prhFileTarget,'file')
   fprintf('Unable to find PRH file %s - check TAG_PATHS and look for file.\n', ...
       prhFileTarget);
   return
end

% load the variables from the file
load(prhFileTarget, 'O', 'K', 'P', 'CUE','TAG','TRIAL');

