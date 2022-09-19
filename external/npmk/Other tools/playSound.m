% Plays a predefined sound. Possible options are:
%
%   failure
%   abort
%   gocue
%   alert
%   warning
%   error
%   done
%   list     - will list names for all the available sounds
%
% Example: playSound('alert')
%
%   Kian Torab
%   kian.torab@utah.edu
%   Department of Bioengineering
%   University of Utah
%   Version 1.2.0 - February 26, 2010

function playSound(soundType)

availableSounds = {'failure',...
                   'abort',...
                   'gocue',...
                   'alert',...
                   'warning',...
                   'done',...
                   'error'};

if ~exist('soundType', 'var')
    disp('Please specify a sound type.');
    return;
end

Freq1   = 0.1;
Freq2   = 0.05;
Freq3   = 0.175;
Freq4   = 0.5;
AmpFrac = 0.5;
FS      = 22050;
%% Failure Sound
FailureSound = [AmpFrac*sin(0:Freq1:Freq1*4000),...
                AmpFrac*sin(0:Freq2:Freq2*10000)]; 
%% Abort Sound
AbortSound   =  AmpFrac*sin(0:Freq2:Freq2*10000);
%% GoCue Sound
GoSound      =  AmpFrac*sin(0:Freq3:Freq3*1000); 
%% Alert Sound
AlertSound   =  AmpFrac*sin(0:Freq4:Freq4*100);
%% Warning Sound
WarningSound =  AmpFrac*sin(0:Freq4:Freq3*10000);
%% Error Sound
TempSound1 = AmpFrac * sin(0:Freq2*1:Freq2*1*6000);
TempSound2 = AmpFrac * sin(0:Freq2*2:Freq2*2*6000);
TempSound3 = AmpFrac * sin(0:Freq2*3:Freq2*3*6000);
ErrorSound = horzcat(TempSound1, TempSound2, TempSound3);
%% Done Sound
quarterBeat = 0.018*6000;
baseFreq = 0.015;
TempSound1 = AmpFrac * sin(0:baseFreq+0.070*1:quarterBeat*4);
TempSound2 = AmpFrac * sin(0:baseFreq+0.050*1:quarterBeat*2);
TempSound3 = AmpFrac * sin(0:baseFreq+0.050*1:quarterBeat*2);
TempSound4 = AmpFrac * sin(0:baseFreq+0.057*1:quarterBeat*4);
TempSound5 = AmpFrac * sin(0:baseFreq+0.050*1:quarterBeat*8);
TempSound6 = AmpFrac * sin(0:baseFreq+0.065*1:quarterBeat*4);
TempSound7 = AmpFrac * sin(0:baseFreq+0.070*1:quarterBeat*8);
DoneSound = horzcat(TempSound1, TempSound2, TempSound3,...
                    TempSound4, TempSound5, TempSound6,...
                    TempSound7);

switch lower(soundType)
    case 'failure'
        playSound = FailureSound;
    case 'abort'
        playSound = AbortSound;
    case 'gocue'
        playSound = GoSound;
    case 'alert'
        playSound = AlertSound;
    case 'warning'
        playSound = WarningSound;
    case 'done'
        playSound = DoneSound;
    case 'error'
        playSound = ErrorSound;
    case 'list'
        disp(availableSounds);
        return;
    otherwise
        disp('Sound type is not valid. Use ''help playSound'' for more information.');
        return;
end

wavplay(playSound, FS);