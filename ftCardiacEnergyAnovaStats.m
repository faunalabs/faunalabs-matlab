%% ftChapterThreeStats.m



%%

figure('Color', white);
[pAovKE_epoch_id, tblAovKE_epoch_id, statsAovKE_epoch_id, ...
    termsAovKE_epoch_id] = ...
    anovan(CC_HR_DATA.iTotalKE, {CC_HR_DATA.epoch CC_HR_DATA.altId}, ...
    'model','interaction','varnames',{'Epoch','ID'} );

[cAovKE_epoch_id, mAovKE_epoch_id, hAovKE_epoch_id, ...
    gnamesAovKE_epoch_id] = ...
    multcompare(statsAovKE, 'Dimension', [1 2] );

figure('Color', white);
boxplot( CC_HR_DATA.iTotalKE, {CC_HR_DATA.id CC_HR_DATA.epoch}, ...
    'Notch','on');


%%

figure; 

s6a = subplot(311);
plot(dT, dFiltLed3*-1, 'r-'); 
xlabel('Time, local');
ylabel('dFiltLed3');
grid; 

s6b = subplot(312); 
plot(dT, dBpAx, 'b-'); 
xlabel('Time, local');
ylabel('bpAx, m · s^2');
grid; 

s6c = subplot(313);
bar(CC_HR_DATA.Time(analysisRange), CC_HR_DATA.iRotKE_x(analysisRange) );
xlabel('Time, local');
ylabel('iTotalKE, J·s');
linkaxes([s6a s6b s6c],'x');