function res = type(this, value)
% Method for and getting/setting EEG file type
% FORMAT res = type(this, value)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: type.m 1267 2008-03-28 12:12:14Z vladimir $

if nargin == 1
    res = this.type;
else
    switch value
        case 'continuous'
            if ntrials(this)>1
                error('Continuous file can only have one trial');
            end
        case 'single'
        case 'evoked' % Add additional checks here
        case 'grandmean'
        otherwise
            error('Unrecognized type');
    end

    this.type = value;
    res = this;
end