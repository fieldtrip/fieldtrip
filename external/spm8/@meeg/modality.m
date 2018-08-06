function [res, list] = modality(this, scalp, planar)
% Returns data modality (like in SPM5)
% FORMAT [res, list] = modality(this, scalp)
%
% scalp - 1 (default) only look at scalp modalities
%         0  look at all modalities
% planar - 1 distinguish between MEG planar and other MEG 
%          0 (default) do not distinguish
% If more than one modality is found the function returns 'Multimodal'
% in res and a cell array of modalities in list.
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: modality.m 2737 2009-02-12 12:49:45Z vladimir $

% Unlike in SPM5, in SPM8 modality is not well defined and is a property
% of channels rather than the whole file. So this function is only a
% temporary solution to make some pieces of code work.

if nargin == 1
    scalp = 1;
end
if nargin < 3
    planar = 0;
end

list = {};

if ~isempty(strmatch('MEG', chantype(this)))
    if planar
        if ~isempty(strmatch('MEGPLANAR', chantype(this)))
            list = [list {'MEGPLANAR'}];
        end
        if ~isempty(strmatch('MEGGRAD', chantype(this))) ||...
                ~isempty(strmatch('MEGMAG', chantype(this))) ||...
                ~isempty(strmatch('MEG', chantype(this), 'exact'))
            list = [list {'MEG'}];
        end
    else
        list = [list {'MEG'}];
    end
end

if ~isempty(strmatch('EEG', chantype(this), 'exact'))
    list = [list {'EEG'}];
end

if ~isempty(strmatch('LFP', chantype(this), 'exact')) && ~scalp
    list = [list {'LFP'}];
end

switch numel(list)
    case 0
        res = 'Other';
    case 1
        res = list{1};
    otherwise
        res = 'Multimodal';
end