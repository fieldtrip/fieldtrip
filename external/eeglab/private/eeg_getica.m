% eeg_getica() - get ICA component activation. Recompute if necessary.
%
% >> mergelocs = eeg_getica(EEG, comp);
%
% Inputs: 
%     EEG     - EEGLAB dataset structure
%     comp    - component index
%
% Output: 
%     icaact  - ICA component activity
%
% Author: Arnaud Delorme, 2006

% Copyright (C) Arnaud Delorme, CERCO, 2006, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function icaact = eeg_getica(EEG, comp)

  if nargin < 1
    help eeg_getica;
    return;
  end
  if nargin < 2
    comp = 1:size(EEG.icaweights,1);
  end
  
  if ~isempty(EEG.icaact)
    icaact = EEG.icaact(comp,:,:);
  else
    disp('Recomputing ICA activations');
    if isempty(EEG.icachansind)
        EEG.icachansind = 1:EEG.nbchan;
        disp('Channels indices are assumed to be in regular order and arranged accordingly');
    end
    icaact = (EEG.icaweights(comp,:)*EEG.icasphere)*reshape(EEG.data(EEG.icachansind,:,:), length(EEG.icachansind), EEG.trials*EEG.pnts);
    icaact = reshape( icaact, size(icaact,1), EEG.pnts, EEG.trials);
  end
