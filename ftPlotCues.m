function numCuePlots = ftPlotCues(cueTimes, cueTypes, yPos)
%
%    numCuePlots = ftPlotCues(cueTimes, cueTypes, yPos);
%
%    Plot CUE.times and CUE.types at yPos on the y-axis.
%    Assumes you're calling this from a time series plotting code block.
%
%    Dave Haas
%    dave.haas@faunalabs.com
%    24 July 2021

if nargin<1
   help ftPlotCues
   return
end

% named colors
blue = [0 0.4470 0.7410];           % '#0072BD';
red = [0.8500 0.3250 0.0980];       % '#D95319'	
goldenrod = [0.9290 0.6940 0.1250];	% '#EDB120'	
purple = [0.4940 0.1840 0.5560];    % '#7E2F8E'	
green = [0.4660 0.6740 0.1880];     % '#77AC30'	
cyan = [0.3010 0.7450 0.9330];      % '#4DBEEE'	
maroon = [0.6350 0.0780 0.1840];    % '#A2142F'
black = [0 0 0];                    % '#000000'
white = [1 1 1];                    % '#FFFFFF'
grey = [0.6 0.6 0.6];               % '#999999'

numCuePlots = 0;

% chillIndex = find(cueTypes(:) == string('1'));
% preApneaIndex = find(cueTypes(:) == string('2'));
% apneaIndex = find(cueTypes(:) == string('3'));
% postApneaIndex = find(cueTypes(:) == string('4'));
respsIndex = find(cueTypes(:) == string('*'));
numResps = numel(respsIndex);

aoIndex = find(cueTypes(:) == string('ao'));
numAO = numel(aoIndex);

acIndex = find(cueTypes(:) == string('ac'));
numAC = numel(acIndex);

excludeIndex = find(cueTypes(:) == string('x'));
numExcludes = numel(excludeIndex);

ventIndex = find(cueTypes(:) == string('v'));
numVents = numel(ventIndex);

led2Index = find(cueTypes(:) == string('1'));
numLed2 = numel(led2Index);

led3Index = find(cueTypes(:) == string('2'));
numLed3 = numel(led3Index);

bothLedsIndex = find(cueTypes(:) == string('3'));
bothLeds = numel(bothLedsIndex);

noLedsIndex = find(cueTypes(:) == string('z'));
noLeds = numel(noLedsIndex);

probOptVentsIndex = find(cueTypes(:) == string('P'));
probOptVents = numel(probOptVentsIndex);

% if (~isempty(chillIndex))
%     plotChillCue = true;
% else 
%     plotChillCue = false;
% end
% 
% if (~isempty(preApneaIndex))
%     plotPreApneaCue = true;
% else 
%     plotPreApneaCue = false;
% end
% 
% if (~isempty(apneaIndex))
%     plotApneaCue = true;
% else 
%     plotApneaCue = false;
% end
% 
% if (~isempty(postApneaIndex))
%     plotPostApneaCue = true;
% else 
%     plotPostApneaCue = false;
% end

if (~isempty(respsIndex))
    plotRespsCue = true;
else 
    plotRespsCue = false;
end

if (~isempty(aoIndex))
    plotAO = true;
else
    plotAO = false;
end

if (~isempty(acIndex))
    plotAC = true;
else
    plotAC = false;
end

if (~isempty(excludeIndex))
    plotExcludes = true;
else
    plotExcludes = false;
end

if (~isempty(ventIndex))
    plotVents = true;
else
    plotVents = false;
end

if (~isempty(led2Index))
    plotLed2 = true;
else
    plotLed2 = false;
end

if (~isempty(led3Index))
    plotLed3 = true;
else
    plotLed3 = false;
end

if (~isempty(bothLedsIndex))
    plotBothLeds = true;
else
    plotBothLeds = false;
end

if (~isempty(noLedsIndex))
    plotNoLeds = true;
else
    plotNoLeds = false;
end

if (~isempty(probOptVentsIndex))
    plotProbOptVents = true;
else
    plotProbOptVents = false;
end

% if (plotChillCue)
%     plot(cueTimes(chillIndex), yPos, 'Marker', 'd', ...
%         'MarkerFaceColor', blue);
%     xline(cueTimes(chillIndex), 'Color', blue, ...
%         'LineStyle', ':', 'LineWidth', 1.2);
%     numCuePlots = numCuePlots + 1;
% end
% 
% if (plotPreApneaCue)
%     plot(cueTimes(preApneaIndex), yPos, 'Marker', '^', ...
%         'MarkerFaceColor', green);
%     xline(cueTimes(preApneaIndex), 'Color', green, ...
%         'LineStyle', ':', 'LineWidth', 1.2);
%     numCuePlots = numCuePlots + 1;
% end
% 
% if (plotApneaCue)
%     plot(cueTimes(apneaIndex), yPos, 'Marker', 'v', ...
%         'MarkerFaceColor', red);
%     xline(cueTimes(apneaIndex), 'Color', red, ...
%         'LineStyle', ':', 'LineWidth', 1.2);
%     numCuePlots = numCuePlots + 1;
% end
% 
% if (plotPostApneaCue)
%     plot(cueTimes(postApneaIndex), yPos, 'Marker', '^', ...
%         'MarkerFaceColor', goldenrod);
%     xline(cueTimes(postApneaIndex), 'Color', goldenrod, ...
%         'LineStyle', ':', 'LineWidth', 1.2);
%     numCuePlots = numCuePlots + 1;
% end

if (plotRespsCue)
    for i = 1:numResps
        plot(cueTimes(respsIndex(i)), yPos, ...
           'Marker', '*', 'MarkerEdgeColor', maroon);
       numCuePlots = numCuePlots + 1;
    end
end

if (plotAO)
   for i = 1:numAO
        plot(cueTimes(aoIndex(i)), yPos, ...
            'Marker', 'd', 'MarkerFaceColor', red, ...
            'MarkerEdgeColor', black);
        xline(cueTimes(aoIndex), 'Color', red, ...
            'LineStyle', ':', 'LineWidth', 1.5);
   end
end

if (plotAC)
   for i = 1:numAC
        plot(cueTimes(acIndex(i)), yPos, ...
            'Marker', 'd', 'MarkerFaceColor', cyan, ...
            'MarkerEdgeColor', black);
        xline(cueTimes(acIndex), 'Color', cyan, ...
            'LineStyle', ':', 'LineWidth', 1.5);
   end
end

if (plotExcludes)
   for i = 1:numExcludes
        plot(cueTimes(excludeIndex(i)), yPos, ...
            'Marker', 'x', 'MarkerFaceColor', purple, ...
            'MarkerEdgeColor', black);
        xline(cueTimes(excludeIndex), 'Color', purple, ...
            'LineStyle', ':', 'LineWidth', 1.5);
   end
end

if (plotVents)
   for i = 1:numVents
        plot(cueTimes(ventIndex(i)), yPos, ...
            'Marker', 'v', 'MarkerFaceColor', green, ...
            'MarkerEdgeColor', black);
        xline(cueTimes(ventIndex), 'Color', green, ...
            'LineStyle', '--', 'LineWidth', 1.5);
   end
end

if (plotLed2)
   for i = 1:numLed2
        plot(cueTimes(led2Index(i)), yPos, ...
            'Marker', 'o', 'MarkerFaceColor', blue, ...
            'MarkerEdgeColor', black);
        xline(cueTimes(led2Index), 'Color', blue, ...
            'LineStyle', '--', 'LineWidth', 1.5);
   end
end

if (plotLed3)
   for i = 1:numLed3
        plot(cueTimes(led3Index(i)), yPos, ...
            'Marker', 'o', 'MarkerFaceColor', red, ...
            'MarkerEdgeColor', black);
        xline(cueTimes(led3Index), 'Color', red, ...
            'LineStyle', '--', 'LineWidth', 1.5);
   end
end

if (plotBothLeds)
   for i = 1:bothLeds
        plot(cueTimes(bothLedsIndex(i)), yPos, ...
            'Marker', 'p', 'MarkerFaceColor', purple, ...
            'MarkerEdgeColor', black, 'MarkerSize', 10);
        xline(cueTimes(bothLedsIndex), 'Color', purple, ...
            'LineStyle', '--', 'LineWidth', 1.5);
   end
end

if (plotNoLeds)
   for i = 1:noLeds
        plot(cueTimes(noLedsIndex(i)), yPos, ...
            'Marker', 's', 'MarkerFaceColor', grey, ...
            'MarkerEdgeColor', black);
        xline(cueTimes(noLedsIndex), 'Color', grey, ...
            'LineStyle', '--', 'LineWidth', 1.5);
   end
end

if (plotProbOptVents)
   for i = 1:probOptVents
        plot(cueTimes(probOptVentsIndex(i)), yPos, ...
            'Marker', '<', 'MarkerFaceColor', cyan, ...
            'MarkerEdgeColor', black);
        xline(cueTimes(probOptVentsIndex), 'Color', cyan, ...
            'LineStyle', '-.', 'LineWidth', 1.5);
   end
end