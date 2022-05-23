function [OPTICS, KINEMATICS, PRESSURE, MAG, POWER, CUE, TAG, TRIAL] = ...
    ftLoadRaw(tag)
%
%    [OPTICS,KINEMATICS,PRESSURE,MAG,POWER,CUE,TAG,TRIAL] = ftLoadRaw(tag);
%
%    Load the post-processed raw sensor data file associated with
%    a tag deployment.
%
%    Dave Haas
%    dave.haas@faunalabs.com
%    21 July 2021

if nargin<1
   help ftLoadRaw
   return
end

OPTICS=[];KINEMATICS=[];PRESSURE=[];MAG=[];POWER=[];CUE=[];TAG=[];TRIAL=[];

global TAG_PATHS;

if isempty(tag)
    fprintf('Specify a valid tag trial, e.g.: tt21_142a and try again.\n');
    return;
end

rawFileTarget = sprintf('%s/%sraw.mat', TAG_PATHS.RAW, tag);

% check if the file exists
if ~exist(rawFileTarget,'file')
   fprintf('Unable to find raw file %s - check TAG_PATHS and look for file.\n', ...
       rawFileTarget);
   return
end

% load the variables from the file
load(rawFileTarget, 'OPTICS', 'KINEMATICS', 'PRESSURE', 'MAG',...
    'POWER','CUE','TAG','TRIAL');

