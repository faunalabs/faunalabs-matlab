%% ftDetrend.m

%   Written by Dave Haas
%   10 July 2021

%   This is an attempt to detrend non-linearities in optical signals

clc;

%%

opticsTime_s = OPTICS.time_s;

figure;
plot(OPTICS.time_s, OPTICS.led2, 'b-', 'LineWidth', 2)
grid on;
erodedSignal = imerode(OPTICS.led2, ones(51, 1));
hold on;
plot(OPTICS.time_s, erodedSignal, 'r-', 'LineWidth', 1);

X = OPTICS.time_s;
Y = erodedSignal;
% Convert X and Y into a table, which is the form fitnlm() likes the input data to be in.
tbl = table(X, Y);

% Define the model as Y = a + exp(-b*x)
% Note how this "x" of modelfun is related to big X and big Y.
% x((:, 1) is actually X and x(:, 2) is actually Y - the first and second columns of the table.

x = zeros(1,length(X));
%modelfun = @(b,x) b(1) + b(2) * exp(-b(3)*x(:, 1));
modelfun = @(b,x) b(1) + b(2) * exp(-b(3)*x(1,:));
beta0 = [10000, 300, 1]; % Guess values to start with.  Just make your best guess.

% Now the next line is where the actual model computation is done.
mdl = fitnlm(tbl, modelfun, beta0);

% Now the model creation is done and the coefficients have been determined.
% YAY!!!!
% Extract the coefficient values from the the model object.
% The actual coefficients are in the "Estimate" column of the "Coefficients" table that's part of the mode.
coefficients = mdl.Coefficients{:, 'Estimate'}
% Create smoothed/regressed data using the model:
yFitted = coefficients(1) + coefficients(2) * exp(-coefficients(3)*X);
% Now we're done and we can plot the smooth model as a red line going through the noisy blue markers.
hold on;
plot(X, yFitted, 'k-', 'LineWidth', 3);
grid on;
title('Exponential Regression with fitnlm()', 'FontSize', fontSize);
xlabel('X', 'FontSize', fontSize);
ylabel('Y', 'FontSize', fontSize);
legendHandle = legend('Original Signal', 'Baseline', 'Fitted Y', 'Location', 'north');
legendHandle.FontSize = 25;
% Set up figure properties:
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');
% Give a name to the title bar.
set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off')