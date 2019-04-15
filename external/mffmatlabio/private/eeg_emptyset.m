% eeg_emptyset() - Initialize an EEG dataset structure with default values.
%
% Usage:
%   >> EEG = eeg_emptyset();
%
% Outputs:
%   EEG    - empty dataset structure with default values.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% 01-25-02 reformated help & license -ad 

function EEG = eeg_emptyset();

EEG.setname     = '';
EEG.filename    = '';
EEG.filepath    = '';
EEG.subject     = '';
EEG.group       = '';
EEG.condition   = '';
EEG.session     = [];
EEG.comments    = '';
EEG.nbchan      = 0;
EEG.trials      = 0;
EEG.pnts        = 0;
EEG.srate       = 1;
EEG.xmin        = 0;
EEG.xmax        = 0;
EEG.times       = [];
EEG.data        = [];
EEG.icaact      = [];
EEG.icawinv     = [];
EEG.icasphere   = [];
EEG.icaweights  = [];
EEG.icachansind = [];
EEG.chanlocs    = [];
EEG.urchanlocs  = [];
EEG.chaninfo    = [];
EEG.ref         = [];
EEG.event       = [];
EEG.urevent     = [];
EEG.eventdescription = {};
EEG.epoch       = [];
EEG.epochdescription = {};
EEG.reject      = [];
EEG.stats       = [];
EEG.specdata    = [];
EEG.specicaact  = [];
EEG.splinefile  = '';
EEG.icasplinefile = '';
EEG.dipfit      = [];
EEG.history     = '';
EEG.saved       = 'no';
EEG.etc         = [];

%EEG.reject.threshold  = [1 0.8 0.85];
%EEG.reject.icareject  = [];
%EEG.reject.compreject = [];
%EEG.reject.gcompreject= [];
%EEG.reject.comptrial  = [];
%EEG.reject.sigreject  = [];
%EEG.reject.elecreject = [];

%EEG.stats.kurta      = [];
%EEG.stats.kurtr      = [];
%EEG.stats.kurtd      = [];		
%EEG.stats.eegentropy = [];
%EEG.stats.eegkurt    = [];
%EEG.stats.eegkurtg   = [];
%EEG.stats.entropy    = [];
%EEG.stats.kurtc      = [];
%EEG.stats.kurtt      = [];
%EEG.stats.entropyc   = [];

return;
