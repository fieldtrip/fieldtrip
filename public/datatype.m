function [type, dimord] = datatype(data, desired)

% DATATYPE determines the type of data represented in a FieldTrip data
% structure and returns a string with raw, freq, timelock source, comp,
% spike, source, volume, dip.
%
% Use as
%   [type, dimord] = datatype(data)
%   [type, dimord] = datatype(data, desired)
%
% See also CHANTYPE, FILETYPE, SENSTYPE, VOLTYPE

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: datatype.m,v $
% Revision 1.8  2009/08/17 08:41:00  jansch
% added undocumented and experimental datatypes mvar and freqmvar
%
% Revision 1.7  2009/06/15 12:55:18  roboos
% also recognize cohspctrm as type=freq
%
% Revision 1.6  2009/03/25 20:52:07  jansch
% changed the conditional order to ensure correct behaviour for a comp-struct
%
% Revision 1.5  2009/01/28 15:00:45  roboos
% fixed detection of timelock for covariance
%
% Revision 1.4  2009/01/28 14:08:30  roboos
% added )
%
% Revision 1.3  2009/01/28 14:08:08  roboos
% return 'unknown'
% detect timelock in case of only data.trial or data.cov
%
% Revision 1.2  2008/12/19 09:12:54  roboos
% added support for desired type, returning boolean
%
% Revision 1.1  2008/12/18 15:49:55  roboos
% new implementation
%

% determine the type of input data, this can be raw, freq, timelock, comp, spike, source, volume, dip
israw      = isfield(data, 'label') && isfield(data, 'time') && isa(data.time, 'cell') && isfield(data, 'trial') && isa(data.trial, 'cell');
isfreq     = isfield(data, 'label') && isfield(data, 'freq') && (isfield(data, 'powspctrm') || isfield(data, 'crsspctrm') || isfield(data, 'cohspctrm') || isfield(data, 'fourierspctrm'));
istimelock = isfield(data, 'label') && isfield(data, 'time') && ~isfield(data, 'freq') && ((isfield(data, 'avg') && isnumeric(data.avg)) || (isfield(data, 'trial') && isnumeric(data.trial) || (isfield(data, 'cov') && isnumeric(data.cov))));
iscomp     = isfield(data, 'topo') || isfield(data, 'topolabel');
isspike    = isfield(data, 'label') && isfield(data, 'waveform') && isa(data.waveform, 'cell') && isfield(data, 'timestamp') && isa(data.timestamp, 'cell');
isvolume   = isfield(data, 'transform') && isfield(data, 'dim');
issource   = isfield(data, 'pos');
isdip      = isfield(data, 'dip');
ismvar     = isfield(data, 'dimord') && ~isempty(strfind(data.dimord, 'lag'));
isfreqmvar = isfield(data, 'freq') && isfield(data, 'transfer');

if iscomp
  type = 'comp';  
  %comp should conditionally go before raw, otherwise the returned datatype
  %will be raw
elseif isfreqmvar
  type = 'freqmvar';
  %freqmvar should conditionally go before freq, otherwise the returned datatype
  %will be freq in the case of frequency mvar data
elseif ismvar
  type = 'mvar';
elseif israw
  type = 'raw';
elseif isfreq
  type = 'freq';
elseif istimelock
  type = 'timelock';
elseif isspike
  type = 'spike';
elseif isvolume
  type = 'volume';
elseif issource
  type = 'source';
elseif isdip
  type = 'dip';
else
  type = 'unknown';
end

if nargin>1
  % return a boolean value
  type = strcmp(type, desired);
  return;
end

if nargout>1
  % also return the dimord of the input data
  if isfield(data, 'dimord')
    dimord = data.dimord;
  else
    dimord = 'unknown';
  end
end

