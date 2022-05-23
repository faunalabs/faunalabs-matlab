function     [RES,AXs,hhh] = fAudit(tag,tcue,RES)
%
%     R = FaunaTagAudit(tag,tcue,RES)
%     Audit tool for FaunaTag
%     tag = tag deployment string e.g., 'ft111_20210514_090956-10.csv'
%     tcue is the time in seconds-since-tag-on to start displaying from
%     RES is an optional audit structure to edit or augment

%     Output:
%        R is the audit structure made in the session. While using this
%        tool, use 's' to save this to a file.
%
%     OPERATION
%
%     FaunaTagAudit will display 20 seconds of FaunaTag data at a time
%     Display from top to bottom is, by default:
%     - Bio-Optics
%     - Aw (m/s^2)
%     - ODBA
%     - Depth
%
%     Within this tool, type or click on the display for these functions:
%
%     - click on the graph to get the time cue, depth, time-to-last
%       and frequency of an event. Time-to-last is the elapsed time 
%       between the current click point and the point last clicked. 
%       Results display in the matlab command window.
%
%     - 'f' to go to the next block
%     - 'b' to go to the previous block
%     - 'j' to jump to a specific time cue. Enter the time cue when
%           prompted, then hit enter.
%     - 'a' to select the current segment and add it to the audit.
%           You will be prompted to enter a cueType on the matlab command
%           window. Enter a single word and type return when complete.
%           E.g.: 'breathHold', 'reorient', 'movement', etc.
%     - 'c' to select the currect cursor position and add it to the 
%           audit as a 0-length comment. Great for marking events, e.g.:
%           TAGON, TAGOFF, RESP, etc.        
%     - 'd' to select current segment as a comment, notated as '*cueName'
%     - 'x' to delete the audit entry at the cursor position (helpful
%           if you make a mistake and want to correct it). If there is no
%           audit entry at the cursor, nothing happens. If there is more 
%           than one audit entry overlapping the cursor, one will be 
%           deleted (the first one encountered in the audit structure).
%     - 'z' to zoom in on a selected segment with PSD plots. When zoomed:
%           - 's', 'd', 'c' work as above
%           - any other keys or clicks close the zoomed-in window
%     - 's' to save and re-sort audit cues to the tag's audit file
%     - 'q' to finish auditing. NOTE: this quits without auto-saving, so
%           remember to save any new audits using 'w' before using 'q'!

%     >>> THE FOLLOWING ARE EXPERIMENTAL / MAY NOT BE WORKING YET <<<
%     - type 'n' to jump to the next time cue stored in the tag audit file.
%       If no audit file exists, or if there are no future time cues stored
%       in the audit file, this has no effect.
%     - type 'l' to jump to the last time cue stored in the tag audit file.
%       If no audit file exists, or if there are no previous time cue in
%       the audit file, this has no effect.
%     - type 'r' to calculate received level between cue selection boundary

%     Dave Haas
%     dave.haas@faunalabs.com
%     Written: 3 June 2021, based on tagaudit.m functionality by M. Johnson
%       to include 'n' and 'l' functionality, and to change MJ's 'l'
%       functionality to 'c' (comments in audit file)

numSeconds = 30;    % number of seconds to display by default

if nargin<2 | isempty(RES),
   RES.cue = [] ;
   RES.comment = [] ;
end


if strcmpi(thirdPlot, 'depth')
    % k = loadprh(tag,0,'p','fs') ;
    k = loadprh(tag,'p','fs') ;
else    
    k = loadprh(tag,'p', 'pitch', 'roll', 'head', 'A', 'Aw', 'fs') ;
end

if k==0,
   fprintf('Unable to find a PRH file - continuing without\n') ;
   p = [] ; fs = [] ;
end

% create 1st derivatives of pitch and heading
if ( (strcmpi(thirdPlot, 'depth')) == 0)
    
    d1Pitch = diff(pitch*180/pi);
    d1Pitch(end+1,1) = 0;                     % pad end w/ zero offset
    d1Roll = diff(roll*90/pi);
    d1Roll(end+1,1) = 0;                      % pad end w/ zero offset
    d1Head = diff(head*180/pi);
    d1Head(end+1,1) = 0;                      % pad end w/ zero offset

    % create ODBA and jerk
    odba    = abs(Aw(1:end,1)) + abs(Aw(1:end,2)) + abs(Aw(1:end,3));
    d1odba  = diff(odba);
    d1odba(end+1,1) = 0; % pad d1odba w/ extra zero; diff removes one element
    jerk  = fs * sqrt ( sum ( abs (diff ( Aw(1:end,1:3).^2 )), 2) );
    jerk(end+1) = 0;   % pad jerk w/ extra zero; diff removes one element

end

% check sampling rate
[x,afs] = tagwavread(tag,tcue,0.01) ;
if SOUND_FH > 0,
   [bs as] = butter(6,SOUND_FH/(afs/2),'high') ;
elseif SOUND_FL > 0,
   [bs as] = butter(6,SOUND_FL/(afs/2)) ;
else
   bs = [] ;
end

% high pass filter for envelope
[bh ah] = cheby1(6,0.5,FH/afs*2,'high') ;
% envelope smoothing filter
pp = 1/TC/afs ;

% angle-of-arrival filter
[baoa aaoa] = butter(4,AOA_FH/(afs/2),'high') ;

current = [0 0] ;
figure(1),clf
if ~isempty(p),
   kb = 1:floor(NS*fs) ;
   AXm = axes('position',[0.11,0.83,0.78,0.10]) ;
   AXc = axes('position',[0.11,0.77,0.78,0.05]) ;
   AXs = axes('position',[0.11,0.45,0.78,0.3]) ;
   AXr = axes('position',[0.11,0.24,0.78,0.17]);
   AXp = axes('position',[0.11,0.10,0.78,0.1]) ;
   

else
   AXm = axes('position',[0.11,0.60,0.78,0.34]) ;
   AXc = axes('position',[0.11,0.52,0.78,0.07]) ;
   AXs = axes('position',[0.11,0.11,0.78,0.38]) ;
end

bc = get(gcf,'Color') ;
set(AXc,'XLim',[0 1],'YLim',[0 1]) ;
set(AXc,'Box','off','XTick',[],'YTick',[],'XColor',bc,'YColor',bc,'Color',bc) ;
cleanh = [] ;

while 1,
   [x,afs] = tagwavread(tag,tcue,NS) ;
   if size(x,2)==1 & nargin>3 & ~isempty(DMON),
      [x,cleanh] = dmoncleanup(x,0.0001,[],cleanh) ;
   end
   if isempty(x), return, end    
   % [B F T] = specgram(x(:,CH),BL,afs,hamming(BL),BL) ;
   [B F T] = specgram(x(:,CH),BL,afs,hamming(BL),round(BL/1.5)) ;
   xx = filter(pp,[1 -(1-pp)],abs(filter(bh,ah,x(:,CH)))) ;

   % plot power
   kk = 1:5:length(xx) ;
   axes(AXm), plot(tcue+kk/afs,xx(kk),'k') ; grid
  
   % Experimental - code to plot upper three quartiles in magenta
   cutOff = 0.003;  % any values <= cutOff are first quartile
   moreThan = xx(kk) > cutOff;
   overX    = tcue+kk(moreThan)/afs;
   overY    = xx(kk(moreThan));
   hold on;
   plot(overX, overY,'m:','MarkerSize',2);
   hold off;
   % End experimental code to make upper three quartiles plot in magenta
   
   set(AXm,'XAxisLocation','top') ;
   yl = get(gca,'YLim') ;
   % Let waveform intensity plot vary from page to page
   % yl(2) = min([yl(2) MAXYONCLICKDISPLAY]) ;
   % Or hard-code from 0 to 0.012, instead of using varying display above
   yl(2) = 0.012;
   axis([tcue tcue+NS yl]) ;
   
   % plot audit cue descriptors
   plotRES(AXc,RES,[tcue tcue+NS]) ;

   % plot depth --or--- kinematic data
   
   if ~isempty(p)
      ks = kb + round(tcue*fs) ;
      
      switch thirdPlot
          
          case 'depth'
              % plot depth
              axes(AXp),plot(ks/fs,p(ks)), grid
              ylabel('Depth, m');
          
          case 'pitch'  
              % plot pitch
              axes(AXp),plot(ks/fs,pitch(ks)*180/pi), grid
              ylabel('Pitch, degrees');
              
          case 'roll'
              % plot roll
              axes(AXp),plot(ks/fs,roll(ks)*90/pi), grid
              ylabel('Roll, degrees');
              
          case 'head'
              % plot head
              axes(AXp),plot(ks/fs,head(ks)*180/pi), grid
              ylabel('Head, degrees');              
              
          case 'd1Pitch' 
              % plot first derivative of pitch
              axes(AXp),plot(ks/fs,d1Pitch(ks)), grid
              ylabel('\deltaPitch (\delta°)');
              
          case 'd1Roll' 
              % plot first derivative of roll
              axes(AXp),plot(ks/fs,d1Roll(ks)), grid
              ylabel('\deltaRoll (\delta°)');
              
          case 'd1Head'       
              % plot first derivative of head
              axes(AXp),plot(ks/fs,d1Head(ks)), grid
              ylabel('\deltaHead (\delta°)');
              
          case 'ODBA'  
              % plot ODBA
              axes(AXp),plot(ks/fs,odba(ks)), grid
              ylabel('ODBA')
              
          case 'd1ODBA'  
              % plot first derivative of ODBA
              axes(AXp),plot(ks/fs,d1odba(ks)), grid
              ylabel('\deltaODBA')
              
          case 'jerk'           
              % plot jerk
              axes(AXp),plot(ks/fs,jerk(ks)), grid
              ylabel('Jerk')
              
          otherwise   
              % plot depth
              axes(AXp),plot(ks/fs,p(ks)), grid
              ylabel('Depth, m');
    
      end
      
   	set(gca,'YDir','reverse') ;
      axis([tcue tcue+max(T) get(gca,'YLim')]) ;
      xlabel('Time, s')
      
   end
   
   % plot spectrogram
   BB = adjust2Axis(20*log10(abs(B))) ;
   colormap default;
   % colormap jet;
   % Default captures full acoustic frequency range
   axes(AXs), imagesc(tcue+T,F/1000,BB,CLIM/1.5) ;
   % This displays the frequency range from 0 - 10 kHz
   % axes(AXs), imagesc(tcue+T,F(1:855)/1000,BB(1:42,:),CLIM/1.5) ;
   axis xy, grid ;
   if ~isempty(p),
      set(AXs,'XTickLabel',[]) ;
   else
      xlabel('Time, s')
   end
   ylabel('Frequency, kHz');
   % Uncomment following two lines to adjust y-axis (frequency)
   % syl = [0 10000];
   % axis([tcue tcue+NS syl]) ;
   hold on
   hhh = plot([0 0],0.8*afs/2000*[1 1],'k*-') ;    % plot cursor
   hold off

   fprintf('\n\n-------------------------\nTime cue = %g\n', tcue);
   
   % Plot rainbow window here
   
   axes(AXr), [cl,a,m,q]=rainbowAlt(x,afs);
   minAoA = min(real(a));
   maxAoA = max(real(a));
   fprintf('MinAoA: %g | MaxAoA: %g\n', minAoA, maxAoA);
   if (minAoA == maxAoA)
       % too few accepted points, so override minAoA and maxAoA
       fprintf('Too few accepted points. Overriding minAoA and maxAoA.\n');
       minAoA = -90;
       maxAoA = 90;
   end
   axis([0 NS minAoA maxAoA]) ;
   ylabel('AoA (Degrees)');
   
   % get input
   done = 0 ;
   while done == 0,
      axes(AXs) ; pause(0) ;
      [gx gy button] = ginput(1) ;
      if button>='A',
         button = lower(setstr(button)) ;
      end
      if button==3 | button=='q',
         save tagaudit_RECOVER RES;
         disp( ' Ending audit program.')
         disp( ' Save your audit data using saveaudit(tag,R) if needed!')
         return

      elseif button=='s',
         ss = input('Enter audit description:','s') ;
         cc = sort(current) ;
         RES.cue = [RES.cue;[cc(1) diff(cc)]] ;
         RES.stype{size(RES.cue,1)} = ss ;
         save tagaudit_RECOVER RES
         plotRES(AXc,RES,[tcue tcue+NS]) ;
      
      elseif button=='d',
         ss = input('Enter non-focal audit description:','s') ;
         ss = sprintf('*%s', ss);
         cc = sort(current) ;
         RES.cue = [RES.cue;[cc(1) diff(cc)]] ;
         RES.stype{size(RES.cue,1)} = ss ;
         save tagaudit_RECOVER RES
         plotRES(AXc,RES,[tcue tcue+NS]) ;
         
      elseif button=='j',
          jj = input('Enter time cue to jump to...');
          tcue = jj - 0.5;
          done = 1;
          
      elseif button=='v',
          vv = input('Enter new volume level (default was 100)...');
          volume = vv;
          done = 1;
          
      elseif button=='w',                    
         saveaudit(tag,RES);
         fprintf('Wrote audit data to %saud.txt\n', tag)
      
      elseif button=='z',                     % -- ZOOM FEATURE STARTS HERE
          % Zoom in on a specific part of the spectrogram
           zoomed = 0; pause(0);
           figure;
           
           ZBL = 8192;

           if (current(2) > current(1))
                [ZZ,Zafs] = tagwavread(tag,current(1),current(2)-current(1)) ;
           elseif (current(1) > current(2))
                [ZZ,Zafs] = tagwavread(tag,current(2),current(1)-current(2)) ;        
           end
           
           % compute the spectrogram
           %[Bz Fz Tz] = specgram(ZZ(:,CH),ZBL,Zafs,hamming(ZBL),ZBL/2) ;
           [Bz Fz Tz] = specgram(ZZ(:,CH),ZBL,Zafs,hamming(ZBL),round(ZBL/1.5)) ;
           
           % compute the PSD
           nFFT     = ZBL;            % segment lengths for fft calculation
           winSize  = hanning(nFFT);  % set hanning window shape
           %nOverlap = nFFT / 2;       % set 50% overlap between segments
           nOverlap = round(nFFT / 1.5);       % set 75% overlap between segments
           [psd_ZZ1, f_ZZ1] = pwelch(ZZ(:,1), winSize, nOverlap, Zafs);
           [psd_ZZ2, f_ZZ2] = pwelch(ZZ(:,2), winSize, nOverlap, Zafs);
           
           % build the zoom plots
           
           AXzw   = axes('position',[0.11,0.9,0.78,0.08]);
           AXzs   = axes('position',[0.11,0.55,0.78,0.30]) ;
           AXzr   = axes('position',[0.11,0.33,0.78,0.17]) ;
           AXpsd  = axes('position',[0.11,0.07,0.78,0.21]) ;
           
           zoomDuration = length(ZZ)/Zafs;
           
           % Plot waveform
           
           waveTime = [1/Zafs:1/Zafs:length(ZZ)/Zafs];
           axes(AXzw), plot(waveTime, ZZ(:,1),'b-');
           hold on; plot(waveTime, ZZ(:,2),'g-');
           xlim([0 zoomDuration]);
           axis xy, grid minor;
           
           % Plot spectrogram
           ccz = sort(current)-tcue ;
           BBz = adjust2Axis(20*log10(abs(Bz))) ;
           % colormap jet;
           colormap default;
           axes(AXzs), imagesc(ccz,Fz/1000,BBz,CLIM) ;
           % This displays the frequency range from 0 - 10 kHz
           % axes(AXzs), imagesc(ccz,Fz(1:428)/1000,BBz(1:41,:),CLIM) ;
           axis xy, grid ;
           %if ~isempty(p),
           %   set(AXzs,'XTickLabel',[]) ;
           %else
           %   xlabel('Time, s')
           %end
           ylabel('Frequency, kHz');
           
           % Plot rainbow window here
           axes(AXzr), [clz,az,mz,qz]=rainbowAlt(ZZ,Zafs);
           minAoAz = min(real(az));
           maxAoAz = max(real(az));
           fprintf('MinAoA: %g | MaxAoA: %g\n', minAoAz, maxAoAz);
           if (minAoAz == maxAoAz)
               % too few accepted points, so override minAoA and maxAoA
               fprintf('Too few accepted points. Overriding minAoA and maxAoA.\n');
               minAoAz = -90;
               maxAoAz = 90;
           end
           maxRainZ = ccz(2) - ccz(1);
           axis([0 maxRainZ minAoAz-10 maxAoAz+10]) ;
           xlabel('Zoom Duration (s)');
           ylabel('AoA (Degrees)');

           % plot power spectral density
           axes(AXpsd), plot(0:Zafs/2, 20*log10(psd_ZZ1), 'b-'); 
           hold on; plot(0:Zafs/2, 20*log10(psd_ZZ2), 'g-');
           xlim([0 Zafs/2]);
           axis xy, grid minor; 
           
           xlabel('Frequency (kHz)');
           ylabel('Power/Frequency (db/Hz)');
           
           while zoomed == 0,               % --- handle zoom window input

               [Zx Zy Zbutton] = ginput(1) ;

               if Zbutton=='e',
                   % export zoomed cue as .wav file to working directory
                   fileBase = tag;
                   fileCue1 = sprintf('_%g',current(1));
                   fileCue2 = sprintf('_%g',current(2));
                   fileExt  = '.wav';
                   fileName = [ fileBase, fileCue1, fileCue2, fileExt ];
                   audiowrite(fileName, ZZ, Zafs);
                   zoomed = 0;
                   
               elseif Zbutton=='p',
                    % play zoomed cue sound
                    chk = min(size(ZZ,2),2) ;
                    if ~isempty(bs),
                       ZZf = filter(bs,as,ZZ(:,1:chk)) ;
                       if ispc
                            sound(volume*ZZf,Zafs/SOUND_DF,16) ;
                       elseif ismac
                            audiowrite('temp.wav', volume*ZZf, Zafs);
                            !afplay temp.wav;
                       end
                    else
                       if ispc
                           sound(volume*ZZ(:,1:chk),Zafs/SOUND_DF,16) ;
                       elseif ismac
                           audiowrite('temp.wav', ZZ(:,1:chk), Zafs);
                           !afplay temp.wav;
                       end
                    end
                    zoomed = 0;

               elseif Zbutton=='s'
                    close;
                    zoomed = 1;
                    ss = input('Enter audit description:','s') ;
                    cc = sort(current) ;
                    RES.cue = [RES.cue;[cc(1) diff(cc)]] ;
                    RES.stype{size(RES.cue,1)} = ss ;
                    save tagaudit_RECOVER RES
                    plotRES(AXc,RES,[tcue tcue+NS]) ;
                   
               elseif Zbutton=='d'     
                    close;
                    zoomed = 1;
                    ss = input('Enter non-focal audit description:','s') ;
                    ss = sprintf('*%s', ss);
                    cc = sort(current) ;
                    RES.cue = [RES.cue;[cc(1) diff(cc)]] ;
                    RES.stype{size(RES.cue,1)} = ss ;
                    save tagaudit_RECOVER RES
                    plotRES(AXc,RES,[tcue tcue+NS]) ;
                                       
               elseif Zbutton=='q',
                   % end zoom, return to runAuditRainbow main routine
                   close;
                   zoomed = 1;
                   
               else
                   % end zoom, return to runAuditRainbow main routine
                   close;
                   zoomed = 1;
               end
           end                                  % -- ZOOM FEATURE ENDS HERE

      elseif button=='a',
         if size(x,2)>1,
            cc = sort(current)-tcue ;
            kcc = round(afs*cc(1)):round(afs*cc(2)) ;
            xf = filter(baoa,aaoa,x(kcc,:)) ;
            [aa,qq] = xc_tdoa(xf(:,1),xf(:,2)) ;
            fprintf(' Angle of arrival %3.1f, quality %1.2f\n',asin(aa*AOA_SCF/afs)*180/pi,qq) ;
         end

      elseif button=='c',
         ss = input(' Enter comment... ','s') ;
         RES.cue = [RES.cue;[gx 0]] ;
         RES.stype{size(RES.cue,1)} = ss ;
         save tagaudit_RECOVER RES
         plotRES(AXc,RES,[tcue tcue+NS]) ;
         
      elseif button=='x',
         kres = min(find(gx>=RES.cue(:,1)-0.1 & gx<sum(RES.cue')'+0.1)) ;
         if ~isempty(kres),
            kkeep = setxor(1:size(RES.cue,1),kres) ;
            RES.cue = RES.cue(kkeep,:) ;
            RES.stype = {RES.stype{kkeep}} ;
            plotRES(AXc,RES,[tcue tcue+NS]) ;
         else
            fprintf(' No saved cue at cursor\n') ;
         end

      elseif button=='f',
            tcue = tcue+floor(NS)-0.5 ;
            done = 1 ;

      elseif button=='b',
            tcue = max([0 tcue-NS+0.5]) ;
            done = 1 ;

      elseif button=='l',
            % set tcue = last time cue stored in tag's audit file
            done = 1;
            
      elseif button=='n',
            % set tcue = next time cue stored in tag's audit file
            done = 1;
            
      elseif button=='r',
            plotRES(AXc,RES,[tcue tcue+NS]) ;
            
      elseif button=='p',
            chk = min(size(x,2),2) ;
            if ~isempty(bs),
               xf = filter(bs,as,x(:,1:chk)) ;
               if ispc
                    sound(volume*xf,afs/SOUND_DF,16) ;
               elseif ismac
                    audiowrite('temp.wav', volume*xf, afs);
                    !afplay temp.wav;
               end
            else
               if ispc
                   sound(volume*x(:,1:chk),afs/SOUND_DF,16) ;
               elseif ismac
                   audiowrite('temp.wav', x(:,1:chk), afs);
                   !afplay temp.wav;
               end
            end

      elseif button==1,
         if gy<0 | gx<tcue | gx>tcue+NS
            fprintf('Invalid click: commands are f b s l p x q\n')

         else
            current = [current(2) gx] ;
            set(hhh,'XData',current) ;
            if ~isempty(p),
               fprintf(' -> %6.1f\t\tdiff to last = %6.1f\t\tp = %6.1f\t\tfreq. = %4.2f kHz\n', ...
                 gx,diff(current),p(round(gx*fs)),gy) ;
			   else
               fprintf(' -> %6.1f\t\tdiff to last = %6.1f\t\tfreq. = %4.2f kHz\n', ...
                 gx,diff(current),gy) ;
	         end
         end
      end
   end
end


function plotRES(AXc,RES,XLIMS) ;
      
axes(AXc)
if ~isempty(RES.cue),
   kk = find(sum(RES.cue')'>XLIMS(1) & RES.cue(:,1)<=XLIMS(2)) ;
   if ~isempty(kk),
      plot([RES.cue(kk,1) sum(RES.cue(kk,:)')']',0.2*ones(2,length(kk)),'k*-') ;
      for k=kk',
         text(max([XLIMS(1) RES.cue(k,1)+0.1]),0.6,RES.stype{k},'FontSize',10) ;
      end
   else
      plot(0,0,'k*-') ;
   end
else
   plot(0,0,'k*-') ;
end

set(AXc,'XLim',XLIMS,'YLim',[0 1]) ;
bc = get(gcf,'Color') ;
set(AXc,'Box','off','XTick',[],'YTick',[],'XColor',bc,'YColor',bc,'Color',bc) ;
return

