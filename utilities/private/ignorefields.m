function ignore = ignorefields(purpose)

% IGNOREFIELDS returns a list of fields that can be present in the cfg structure that
% should be ignored at various places in the code, e.g. for provenance, history,
% size-checking, etc.

switch purpose

  case 'appendtimelock'
    ignore = {
      'cfg'
      'dimord'
      'elec'
      'fsample'
      'grad'
      'label'
      'negclusters'
      'opto'
      'posclusters'
      'sampleinfo' % this is dealt with explicitly
      'time'
      'topo'
      'topodimord'
      'topolabel'
      'trialinfo'  % this is dealt with explicitly
      'unmixing'
      'unmixingdimord'
      };

  case 'appendfreq'
    ignore = {
      'cfg'
      'cumsumcnt'  % this is dealt with explicitly
      'cumtapcnt'  % this is dealt with explicitly
      'dimord'
      'elec'
      'freq'
      'grad'
      'label'
      'negclusters'
      'opto'
      'posclusters'
      'sampleinfo' % this is dealt with explicitly
      'time'
      'trialinfo'  % this is dealt with explicitly
      };

  case 'deface'
    ignore = {
      % some fields should be dealt with explicitly
      'dim'
      'hex'
      'pos'
      'tet'
      'transform'
      'tri'
      % some fields are irrelevant
      'cfg'
      'coordsys'
      'fid'
      'unit'
      };

  case 'pipeline'
    ignore = {
      % some fields that are always allowed to be present in the configuration
      'cfg'
      'inside'
      'leadfield'
      'previous'
      };

  case 'allowed'
    ignore = {
      % some fields that are always allowed to be present in the configuration
      'callinfo'
      'checkconfig'
      'checkpath'
      'checksize'
      'debug'
      'notification'
      'outputfilepresent'
      'previous'
      'progress'
      'showcallinfo'
      'spmversion'
      'toolbox'
      'trackcallinfo'
      'trackconfig'
      'trackdatainfo'
      'trackusage'
      'version'
      'warning'
      };

  case {'rollback'}
    ignore = {
      % these should not be updated in rollback_provenance
      'callinfo'
      'checkconfig'
      'checksize'
      'debug'
      'notification'
      'previous'
      'showcallinfo'
      'trackcallinfo'
      'trackconfig'
      'trackdatainfo'
      'trackusage'
      'version'
      'warning'
      };

  case {'provenance', 'history'}
    ignore = {
      % these should not be included in the provenance or history
      'checkconfig'
      'checksize'
      'debug'
      'notification'
      'reproducescript'
      'showcallinfo'
      'trackcallinfo'
      'trackconfig'
      'trackdatainfo'
      'trackusage'
      'warning'
      };

  case {'reproducescript'}
    ignore = {
      % these should not be included in the output script
      'callinfo'
      'checkconfig'
      'checkpath'
      'checksize'
      'debug'
      'notification'
      'outputfilepresent'
      'progress'
      'reproducescript'
      'showcallinfo'
      'toolbox'
      'trackcallinfo'
      'trackconfig'
      'trackdatainfo'
      'trackusage'
      'version'
      'warning'
      };

  case 'trackconfig'
    ignore = {
      % these fields from the user should be ignored
      'artfctdef'
      'artifact'
      'checksize'
      'event'
      'trl'
      'trlold'
      % these fields are for internal usage only
      'callinfo'
      'checkconfig'
      'checksize'
      'debug'
      'hastoolbox'
      'notification'
      'previous'
      'showcallinfo'
      'trackcallinfo'
      'trackconfig'
      'trackdatainfo'
      'trackusage'
      'version'
      'warning'
      };

  case 'checksize'
    ignore = {
      % the size of these fields should not be checked
      'artfctdef'
      'artifact'
      'checksize'
      'event'
      'hastoolbox'
      'previous'
      'trl'
      'trlold'
      };

  case 'makessense'
    ignore = {
      % these fields should not be used to check whether the trialinfo and sampleinfo make sense
      'cfg'
      'dimord'
      'elec'
      'freq'
      'fsample'
      'grad'
      'hdr'
      'label'
      'negclusters'
      'opto'
      'posclusters'
      'sampleinfo'
      'time'
      'trialinfo'
      };

  case 'html'
    ignore = {
      % when generating a html-formatted pipeline, ignore data-like fields and fields that probably were not added by the user himself
      'checkconfig'
      'checksize'
      'debug'
      'event'
      'headmodel'
      'outputfilepresent'
      'previous'
      'progress'
      'showcallinfo'
      'sourcemodel'
      'trackcallinfo'
      'trackconfig'
      'trackdatainfo'
      'trackusage'
      'warning'
      };

  case 'selectdata'
    ignore = {
      % these fields do not contain data and should be excluded
      'cfg'
      'coordsys'
      'dim'
      'elec'
      'fsample'
      'fsampleorig'
      'grad'
      'hdr'
      'negclusters'
      'opto'
      'posclusters'
      'topolabel'
      'transform'
      'unit'
      };

  case 'recursesize'
    ignore = {
      % these fields should not recursively be checked on their size
      'elec'
      'event'
      'grad'
      'grid'
      'headmodel'
      'headshape'
      'layout'
      'matrix'
      'mri'
      'neighbours'
      'opto'
      'sourcemodel'
      'vol'
      };

  otherwise
    ft_error('invalid purpose');
end % switch purpose

