function res = chanlabels(this, varargin)
% Method for getting/setting the channel labels
% FORMAT res = chanlabels(this, ind, label)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: chanlabels.m 1373 2008-04-11 14:24:03Z spm $

if nargin == 3
    ind = varargin{1};
    label = varargin{2};
    
    if iscell(label) && length(label)>1
        if ~isempty(ind) && length(ind)~=length(label)
            error('Indices and values do not match');
        end
        
        if length(label)>1
            for i = 1:length(label)
                for j = (i+1):length(label)
                    if strcmp(label{i}, label{j})
                        error('All labels must be different');
                    end
                end
            end
        end
        
    end
end

res = getset(this, 'channels', 'label', varargin{:});
