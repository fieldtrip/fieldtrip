function ignore = ignorefields(purpose)

% IGNOREFIELDS returns a list of fields that can be present in the cfg structure that
% should be ignored at various places in the code, e.g. for provenance, history,
% size-checking, etc.

switch purpose
  
  case 'appendtimelock'
    ignore = {
      'cfg'
      'label'
      'time'
      'dimord'
      'grad'
      'elec'
      'opto'
      'trialinfo'  % this is dealt with explicitly
      'sampleinfo' % this is dealt with explicitly
      };
    
  case 'appendfreq'
    ignore = {
      'cfg'
      'label'
      'time'
      'freq'
      'dimord'
      'grad'
      'elec'
      'opto'
      'trialinfo'  % this is dealt with explicitly
      'sampleinfo' % this is dealt with explicitly
      'cumsumcnt'  % this is dealt with explicitly
      'cumtapcnt'  % this is dealt with explicitly
      };
    
  case 'deface'
    ignore = {
      % some fields should be dealt with explicitly
      'pos'
      'tri'
      'tet'
      'hex'
      'dim'
      'transform'
      % some fields are irrelevant
      'unit'
      'coordsys'
      'fid'
      'cfg'
      };
    
  case 'pipeline'
    ignore = {
      % some fields that are always allowed to be present in the configuration
      'leadfield'
      'inside'
      'cfg'
      'previous'
      };
    
  case 'allowed'
    ignore = {
      % some fields that are always allowed to be present in the configuration
      'postamble'
      'trackconfig'
      'checkconfig'
      'checksize'
      'trackusage'
      'trackdatainfo'
      'trackcallinfo'
      'showcallinfo'
      'callinfo'
      'version'
      'warning'
      'notification'
      'debug'
      'previous'
      'progress'
      'outputfilepresent'
      'toolbox'
      };
    
    
  case {'provenance', 'history'}
    ignore = {
      % these should not be included in the provenance or history
      'postamble'
      'checkconfig'
      'checksize'
      'trackconfig'
      'trackusage'
      'trackdatainfo'
      'trackcallinfo'
      'showcallinfo'
      'warning'
      'notification'
      'debug'
      'progress'
      };
    
    
  case 'trackconfig'
    ignore = {
      % these fields from the user should be ignored
      'checksize'
      'trl'
      'trlold'
      'event'
      'artifact'
      'artfctdef'
      % these fields are for internal usage only
      'postamble'
      'checkconfig'
      'checksize'
      'trackconfig'
      'trackusage'
      'trackdatainfo'
      'trackcallinfo'
      'showcallinfo'
      'callinfo'
      'version'
      'warning'
      'notification'
      'debug'
      'previous'
      };
    
  case 'checksize'
    ignore = {
      % the size of these fields should not be checked
      'checksize'
      'trl'
      'trlold'
      'event'
      'artifact'
      'artfctdef'
      'previous'
      };
    
  otherwise
    ft_error('invalid purpose');
end % switch purpose
