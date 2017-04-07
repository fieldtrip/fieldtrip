function res = sensors(this, type, newsens)
% Sets and gets sensor fields for EEG and MEG
% returns empty matrix if no sensors are defined.
% FORMAT res = sensors(this, type, newsens)
%   type - 'EEG' or 'MEG'
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: sensors.m 4790 2012-07-17 16:00:20Z vladimir $

if nargin<2
    error('Sensor type (EEG or MEG) must be specified');
end    

switch lower(type)
    case 'eeg'
        if nargin < 3
            if isfield(this.sensors, 'eeg')
                res = this.sensors.eeg;
            else
                res = [];
            end
        else
            this.sensors(1).eeg = ft_datatype_sens(newsens);
            res = this;
        end
    case 'meg'
        if nargin < 3
            if isfield(this.sensors, 'meg')
                res = this.sensors.meg;
            else
                res = [];
            end
        else
            this.sensors(1).meg = ft_datatype_sens(newsens);
            res = this;
        end
    otherwise
        error('Unsupported sensor type');
end
