%% ftGrammPlots.m

%   written by Dave Haas, 5-7 August 2021

%%  plot using gramm

cN = width(candidateOdbaRidge);
cX = length(candidateOdbaRidge);
cval = {'bpOdba','hrBpOdba','sqwaveBpOdba','sumImfBpOdba'};

kT_x = 1:1:cX;

clear g;
g(1,1)=gramm('x', kT_s, 'y', candidateOdbaRidge', 'color', cval );
g(1,2)=copy(g(1));
g(2,1)=copy(g(1));
g(2,2)=copy(g(1));

g(1,1).geom_point();
g(1,1).set_title('geom_point()');

g(1,2).geom_line();
g(1,2).set_title('geom_line()');

g(2,1).stat_smooth();
g(2,1).set_point_options('base_size',3);
g(2,1).set_title('stat_smooth()');


g(2,2).stat_summary();
g(2,2).set_title('stat_summary()');

g.set_title('Candidate Heart Rate using ODBA-based WSST Ridging ');

figure('Position',[100 100 800 550]);
g.draw();