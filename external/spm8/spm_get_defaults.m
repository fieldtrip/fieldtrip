function varargout = spm_get_defaults(defstr, varargin)
% Get/set the defaults values associated with an identifier
% FORMAT defaults = spm_get_defaults
% Return the global "defaults" variable defined in spm_defaults.m.
%
% FORMAT defval = spm_get_defaults(defstr)
% Return the defaults value associated with identifier "defstr". 
% Currently, this is a '.' subscript reference into the global  
% "defaults" variable defined in spm_defaults.m.
%
% FORMAT spm_get_defaults(defstr, defval)
% Sets the defaults value associated with identifier "defstr". The new
% defaults value applies immediately to:
% * new modules in batch jobs
% * modules in batch jobs that have not been saved yet
% This value will not be saved for future sessions of SPM. To make
% persistent changes, edit spm_defaults.m.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Volkmar Glauche
% $Id: spm_get_defaults.m 3692 2010-01-21 21:43:31Z guillaume $

global defaults;
if isempty(defaults)
    spm_defaults;
end

if nargin == 0
    varargout{1} = defaults;
    return
end

% construct subscript reference struct from dot delimited tag string
tags = textscan(defstr,'%s', 'delimiter','.');
subs = struct('type','.','subs',tags{1}');

if nargin == 1
    varargout{1} = subsref(defaults, subs);
else
    defaults = subsasgn(defaults, subs, varargin{1});
end
