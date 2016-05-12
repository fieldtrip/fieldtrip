function ignore = ignorefields(varargin)

% IGNOREFIELDS returns a list of fields that can be present in the cfg structure 
% that should be ignored for provenance and history.

% this is a list of fields that should not be included in provenance or history
ignore = {
  'trackconfig'
  'trackusage'
  'trackdatainfo'
  'trackcallinfo'
  'showcallinfo'
  'warning'
  'progress'
};

