function [cfg] = ft_checkconfig(cfg, varargin)

% FT_CHECKCONFIG checks the input cfg of the main FieldTrip functions.
%
% 1: It checks whether the cfg contains all the required options, it gives
% a warning when renamed or deprecated options are used, and it makes sure
% no forbidden options are used. If necessary and possible, this function
% will adjust the cfg to the input requirements. If the input cfg does NOT
% correspond to the requirements, this function gives an elaborate warning
% message.
%
% 2: It controls the relevant cfg options that are being passed on to other
% functions, by putting them into substructures or converting them into the
% required format.
%
% 3: It controls the output cfg (data.cfg) such that it only contains
% relevant and used fields. The size of fields in the output cfg is also
% controlled: fields exceeding a certain maximum size are emptied.
% This part of the functionality is still under construction!
%
% Use as
%   [cfg] = ft_checkconfig(cfg, ...)
%
% The behaviour of checkconfig can be controlled by the following cfg options,
% which can be set as global fieldtrip defaults (see FT_DEFAULTS)
%   cfg.checkconfig = 'pedantic', 'loose' or 'silent' (control the feedback behaviour of checkconfig)
%   cfg.trackconfig = 'cleanup', 'report' or 'off'
%   cfg.checksize   = number in bytes, can be inf (set max size allowed for output cfg fields)
%
% Optional input arguments should be specified as key-value pairs and can include
%   renamed         = {'old',  'new'}        % list the old and new option
%   renamedval      = {'opt',  'old', 'new'} % list option and old and new value
%   allowedval      = {'opt', 'allowed1'...} % list of allowed values for a particular option, anything else will throw an error
%   required        = {'opt1', 'opt2', etc.} % list the required options
%   allowed         = {'opt1', 'opt2', etc.} % list the allowed options, all other options are forbidden
%   forbidden       = {'opt1', 'opt2', etc.} % list the forbidden options, these result in an error
%   deprecated      = {'opt1', 'opt2', etc.} % list the deprecated options
%   unused          = {'opt1', 'opt2', etc.} % list the unused options, these will be removed and a warning is issued
%   createsubcfg    = {'subname', etc.}      % list the names of the subcfg
%   dataset2files   = 'yes', 'no'            % converts dataset into headerfile and datafile
%   index2logical   = 'yes', 'no'            % converts cfg.index or cfg.grid.index into logical representation
%   checksize       = 'yes', 'no'            % remove large fields from the cfg
%   trackconfig     = 'on', 'off'            % start/end config tracking
%
% See also FT_CHECKDATA, FT_DEFAULTS

% Copyright (C) 2007-2014, Robert Oostenveld, Saskia Haegens
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

renamed         = ft_getopt(varargin, 'renamed');
allowed         = ft_getopt(varargin, 'allowed');
required        = ft_getopt(varargin, 'required');
deprecated      = ft_getopt(varargin, 'deprecated');
unused          = ft_getopt(varargin, 'unused');
forbidden       = ft_getopt(varargin, 'forbidden');
renamedval      = ft_getopt(varargin, 'renamedval');
allowedval      = ft_getopt(varargin, 'allowedval');
createsubcfg    = ft_getopt(varargin, 'createsubcfg');
checkfilenames  = ft_getopt(varargin, 'dataset2files');
checkinside     = ft_getopt(varargin, 'index2logical', 'off');
checksize       = ft_getopt(varargin, 'checksize', 'off');
trackconfig     = ft_getopt(varargin, 'trackconfig');

if ~isempty(trackconfig) && strcmp(trackconfig, 'on')
  % infer from the user configuration whether tracking should be enabled
  if isfield(cfg, 'trackconfig') && (strcmp(cfg.trackconfig, 'report') || strcmp(cfg.trackconfig, 'cleanup'))
    trackconfig = 'on'; % turn on configtracking if user requests report/cleanup
  else
    trackconfig = []; % disable configtracking if user doesn't request report/cleanup
  end
end

% these should be cell arrays and not strings
if ischar(required),     required     = {required};      end
if ischar(deprecated),   deprecated   = {deprecated};    end
if ischar(unused),       unused       = {unused};        end
if ischar(forbidden),    forbidden    = {forbidden};     end
if ischar(createsubcfg), createsubcfg = {createsubcfg};  end

if isfield(cfg, 'checkconfig')
  silent   = strcmp(cfg.checkconfig, 'silent');
  loose    = strcmp(cfg.checkconfig, 'loose');
  pedantic = strcmp(cfg.checkconfig, 'pedantic');
else
  silent   = false;
  loose    = true;
  pedantic = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rename old to new options, give warning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(renamed)
  if issubfield(cfg, renamed{1})
    cfg = setsubfield(cfg, renamed{2}, (getsubfield(cfg, renamed{1})));
    cfg = rmsubfield(cfg, renamed{1});
    if silent
      % don't mention it
    elseif loose
      warning('use cfg.%s instead of cfg.%s', renamed{2}, renamed{1});
    elseif pedantic
      error('use cfg.%s instead of cfg.%s', renamed{2}, renamed{1});
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rename old to new value, give warning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(renamedval) && issubfield(cfg, renamedval{1})
  if strcmpi(getsubfield(cfg, renamedval{1}), renamedval{2})
    cfg = setsubfield(cfg, renamedval{1}, renamedval{3});
    if silent
      % don't mention it
    elseif loose
      warning('use cfg.%s=''%s'' instead of cfg.%s=''%s''', renamedval{1}, renamedval{3}, renamedval{1}, renamedval{2});
    elseif pedantic
      error('use cfg.%s=''%s'' instead of cfg.%s=''%s''', renamedval{1}, renamedval{3}, renamedval{1}, renamedval{2});
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for required fields, give error when missing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(required)
  fieldsused = fieldnames(cfg);
  [c, ia, ib] = setxor(required, fieldsused);
  if ~isempty(ia)
    error('The field cfg.%s is required\n', required{ia});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for deprecated fields, give warning when present
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(deprecated)
  fieldsused = fieldnames(cfg);
  if any(ismember(deprecated, fieldsused))
    if silent
      % don't mention it
    elseif loose
      warning('The option cfg.%s is deprecated, support is no longer guaranteed\n', deprecated{ismember(deprecated, fieldsused)});
    elseif pedantic
      error('The option cfg.%s is not longer supported\n', deprecated{ismember(deprecated, fieldsused)});
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for unused fields, give warning when present and remove them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(unused)
  fieldsused = fieldnames(cfg);
  if any(ismember(unused, fieldsused))
    cfg = rmfield(cfg, unused(ismember(unused, fieldsused)));
    if silent
      % don't mention it
    elseif loose
      warning('The field cfg.%s is unused, it will be removed from your configuration\n', unused{ismember(unused, fieldsused)});
    elseif pedantic
      error('The field cfg.%s is unused\n', unused{ismember(unused, fieldsused)});
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for required fields, give error when missing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(allowed)
  % there are some general options that should always be allowed
  allowed = union(allowed, {
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
    'debug'
    'previous'
    'outputfilepresent'
    });
  fieldsused = fieldnames(cfg);
  [c, i] = setdiff(fieldsused, allowed);
  if ~isempty(c)
    error('The field cfg.%s is not allowed\n', c{1});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for forbidden fields, give error when present
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(forbidden)
  fieldsused = fieldnames(cfg);
  if any(ismember(forbidden, fieldsused))
    cfg = rmfield(cfg, forbidden(ismember(forbidden, fieldsused)));
    if silent
      % don't mention it
    elseif loose
      warning('The field cfg.%s is forbidden, it will be removed from your configuration\n', forbidden{ismember(forbidden, fieldsused)});
    elseif pedantic
      error('The field cfg.%s is forbidden\n', forbidden{ismember(forbidden, fieldsused)});
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for allowed values, give error if non-allowed value is specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(allowedval) && isfield(cfg, allowedval{1}) ...
    && ~any(strcmp(cfg.(allowedval{1}), allowedval(2:end)))
  s = ['The only allowed values for cfg.' allowedval{1} ' are: '];
  for k = 2:numel(allowedval)
    s = [s allowedval{k} ', '];
  end
  s = s(1:end-2); % strip last comma
  error(s);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% backward compatibility for the gradiometer and electrode definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(cfg, 'grad') && ~isempty(cfg.grad)
  cfg.grad = ft_datatype_sens(cfg.grad);
end
if isfield(cfg, 'elec')&& ~isempty(cfg.elec)
  cfg.elec = ft_datatype_sens(cfg.elec);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% backward compatibility for old neighbourstructures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(cfg, 'neighbours') && iscell(cfg.neighbours)
  warning('Neighbourstructure is in old format - converting to structure array');
  cfg.neighbours= fixneighbours(cfg.neighbours);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% createsubcfg
%
% This collects the optional arguments for some of the low-level
% functions and puts them in a separate substructure. This function is to
% ensure backward compatibility of end-user scripts, fieldtrip functions
% and documentation that do not use the nested detailled configuration
% but that use a flat configuration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(createsubcfg)
  for j=1:length(createsubcfg)
    subname = createsubcfg{j};
    
    if isfield(cfg, subname)
      % get the options that are already specified in the substructure
      subcfg = getfield(cfg, subname);
    else
      % start with an empty substructure
      subcfg = [];
    end
    
    % add all other relevant options to the substructure
    switch subname
      case 'preproc'
        fieldname = {
          'reref'
          'refchannel'
          'implicitref'
          'detrend'
          'bpfiltdir'
          'bpfilter'
          'bpfiltord'
          'bpfilttype'
          'bpfreq'
          'bsfiltdir'
          'bsfilter'
          'bsfiltord'
          'bsfilttype'
          'bsfreq'
          'demean'
          'baselinewindow'
          'denoise'
          'dftfilter'
          'dftfreq'
          'hpfiltdir'
          'hpfilter'
          'hpfiltord'
          'hpfilttype'
          'hpfreq'
          'lpfiltdir'
          'lpfilter'
          'lpfiltord'
          'lpfilttype'
          'lpfreq'
          'medianfilter'
          'medianfiltord'
          'hilbert'
          'derivative'
          'rectify'
          'boxcar'
          'absdiff'
          };
        
      case 'grid'
        fieldname = {
          'xgrid'
          'ygrid'
          'zgrid'
          'resolution'
          'unit'
          'filter'
          'leadfield'
          'inside'
          'outside'
          'pos'
          'dim'
          'tight'
          };
        
      case 'dics'
        fieldname = {
          'feedback'
          'fixedori'
          'keepcsd'
          'keepfilter'
          'keepmom'
          'keepsubspace'
          'lambda'
          'normalize'
          'normalizeparam'
          'powmethod'
          'projectnoise'
          'reducerank'
          'realfilter'
          'subspace'
          };
        
      case 'eloreta'
        fieldname = {
          'keepfilter'
          'keepmom'
          'lambda'
          'normalize'
          'normalizeparam'
          'reducerank'
          };
        
      case 'lcmv'
        fieldname = {
          'feedback'
          'fixedori'
          'keepcov'
          'keepfilter'
          'keepmom'
          'keepsubspace'
          'lambda'
          'normalize'
          'normalizeparam'
          'powmethod'
          'projectnoise'
          'projectmom'
          'reducerank'
          'subspace'
          };
        
      case 'pcc'
        fieldname = {
          'feedback'
          'keepfilter'
          'keepmom'
          'lambda'
          'normalize'
          'normalizeparam'
          %'powmethod'
          'projectnoise'
          'reducerank'
          'keepcsd'
          'realfilter'
          'fixedori'
          };
        
      case {'rv'}
        fieldname = {
          'feedback'
          'lambda'
          };
        
      case 'mne'
        fieldname = {
          'feedback'
          'lambda'
          'keepfilter'
          'prewhiten'
          'snr'
          'scalesourcecov'
          };
        
      case 'music'
        fieldname = {
          'feedback'
          'numcomponent'
          };
        
      case 'sam'
        fieldname = {
          'meansphereorigin'
          'feedback'
          'lambda'
          'fixedori'
          'reducerank'
          'normalize'
          'normalizeparam'
          };
        
      case 'mvl'
        fieldname = {};
        
      case {'npsf', 'granger'}
        % non-parametric spectral factorization -> csd2transfer
        fieldname = {
          'block'
          'blockindx'
          'channelcmb'
          'numiteration'
          'tol'
          'sfmethod'
          'svd'
          'init'
          'checkconvergence'
          };
        
      otherwise
        error('unexpected name of the subfunction');
        fieldname = {};
        
    end % switch subname
    
    for i=1:length(fieldname)
      if ~isfield(subcfg, fieldname{i}) && isfield(cfg, fieldname{i})
        
        if silent
          % don't mention it
        elseif loose
          warning('The field cfg.%s is deprecated, please specify it as cfg.%s.%s instead of cfg.%s', fieldname{i}, subname, fieldname{i});
        elseif pedantic
          error('The field cfg.%s is not longer supported, use cfg.%s.%s instead\n', fieldname{i}, subname, fieldname{i});
        end
        
        subcfg = setfield(subcfg, fieldname{i}, getfield(cfg, fieldname{i}));  % set it in the subconfiguration
        cfg = rmfield(cfg, fieldname{i});                                      % remove it from the main configuration
      end
    end
    
    % copy the substructure back into the main configuration structure
    cfg = setfield(cfg, subname, subcfg);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checkinside, i.e. index2logical
%
% Converts indexed cfg.inside/outside into logical representation if neccessary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if istrue(checkinside)
  if isfield(cfg, 'inside') && any(cfg.inside>1)
    inside = false(size(cfg.pos,1),1);
    inside(cfg.inside) = true;
    cfg = removefields(cfg, {'inside', 'outside'});
    cfg.inside = inside;
  elseif isfield(cfg, 'grid') && isfield(cfg.grid, 'inside') && any(cfg.grid.inside>1)
    inside = false(size(cfg.grid.pos,1),1);
    inside(cfg.grid.inside) = true;
    cfg.grid = removefields(cfg.grid, {'inside', 'outside'});
    cfg.grid.inside = inside;
  end
end % if checkinside

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checkfilenames, i.e. dataset2files
%
% Converts cfg.dataset into cfg.headerfile and cfg.datafile if neccessary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if istrue(checkfilenames)
  
  % start with empty fields if they are not present
  if ~isfield(cfg, 'dataset')
    cfg.dataset = [];
  end
  if ~isfield(cfg, 'datafile')
    cfg.datafile = [];
  end
  if ~isfield(cfg, 'headerfile')
    cfg.headerfile = [];
  end
  
  if ~isempty(cfg.dataset)
    % the dataset is an abstract concept and might relate to a file, a
    % constellation of fioles or a directory containing multiple files
    
    if isequal(cfg.dataset, 'gui') || isequal(cfg.dataset, 'uigetfile')
      % display a graphical file selection dialog
      [f, p] = uigetfile('*.*', 'Select a data file');
      if isequal(f, 0)
        error('User pressed cancel');
      else
        d = fullfile(p, f);
      end
      cfg.dataset = d;
    elseif strcmp(cfg.dataset, 'uigetdir')
      % display a graphical directory selection dialog
      d = uigetdir('*.*', 'Select a data directory');
      if isequal(d, 0)
        error('User pressed cancel');
      end
      cfg.dataset = d;
    end
    
    % ensure that the headerfile and datafile are defined, which are sometimes different than the name of the dataset
    % this requires correct autodetection of the format of the data set
    [cfg.dataset, cfg.headerfile, cfg.datafile] = dataset2files(cfg.dataset, []);
    
  elseif ~isempty(cfg.datafile) && isempty(cfg.headerfile);
    % assume that the datafile also contains the header information
    cfg.dataset    = cfg.datafile;
    cfg.headerfile = cfg.datafile;
    
  elseif isempty(cfg.datafile) && ~isempty(cfg.headerfile);
    % assume that the headerfile also contains the data
    cfg.dataset  = cfg.headerfile;
    cfg.datafile = cfg.headerfile;
  end
  
  % fill dataformat if unspecified, doing this only once saves time later
  if ~isfield(cfg,'dataformat') || isempty(cfg.dataformat)
    cfg.dataformat = ft_filetype(cfg.datafile);
  end
  
  % fill headerformat if unspecified, doing this only once saves time later
  if ~isfield(cfg,'headerformat') || isempty(cfg.headerformat)
    cfg.headerformat = ft_filetype(cfg.headerfile);
  end
  
  % remove empty fields, otherwise a subsequent check on required fields doesn't make any sense
  if isempty(cfg.dataset),    cfg = rmfield(cfg, 'dataset');    end
  if isempty(cfg.headerfile), cfg = rmfield(cfg, 'headerfile'); end
  if isempty(cfg.datafile),   cfg = rmfield(cfg, 'datafile');   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% configtracking
%
% switch configuration tracking on/off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(trackconfig)
  try
    if strcmp(trackconfig, 'on')
      if isa(cfg, 'struct')
        % turn ON configuration tracking
        cfg = config(cfg);
        % remember that configtracking has been turned on
        cfg = access(cfg, 'set', 'counter', 1);
      elseif isa(cfg, 'config')
        % remember how many times trackconfig has been turned on
        cfg = access(cfg, 'set', 'counter', access(cfg, 'get', 'counter')+1); % count the 'ONs'
      end
    end
    
    if strcmp(trackconfig, 'off') && isa(cfg, 'config')
      % turn OFF configuration tracking, optionally give report and/or cleanup
      cfg = access(cfg, 'set', 'counter', access(cfg, 'get', 'counter')-1); % count(down) the 'OFFs'
      
      if access(cfg, 'get', 'counter')==0
        % only proceed when number of 'ONs' matches number of 'OFFs'
        
        if strcmp(cfg.trackconfig, 'report') || strcmp(cfg.trackconfig, 'cleanup')
          % gather information about the tracked results
          r = access(cfg, 'reference');
          o = access(cfg, 'original');
          
          key = fieldnames(cfg);
          key = key(:)';

          ignorefields = {
             % these fields from the user should never be removed
             'trl'
             'trlold'
             'event'
             'artifact'
             'artfctdef'
             % these fields are for internal usage
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
             'debug'
             'previous'
           };

          skipsel      = match_str(key, ignorefields);
          key(skipsel) = [];
          
          used     = zeros(size(key));
          original = zeros(size(key));
          
          for i=1:length(key)
            used(i)     = (r.(key{i})>0);
            original(i) = (o.(key{i})>0);
          end
          
          if ~silent
            % give report on screen
            fprintf('\nThe following config fields were specified by YOU and were USED\n');
            sel = find(used & original);
            if numel(sel)
              fprintf('  cfg.%s\n', key{sel});
            else
              fprintf('  <none>\n');
            end
            
            fprintf('\nThe following config fields were specified by YOU and were NOT USED\n');
            sel = find(~used & original);
            if numel(sel)
              fprintf('  cfg.%s\n', key{sel});
            else
              fprintf('  <none>\n');
            end
            
            fprintf('\nThe following config fields were set to DEFAULTS and were USED\n');
            sel = find(used & ~original);
            if numel(sel)
              fprintf('  cfg.%s\n', key{sel});
            else
              fprintf('  <none>\n');
            end
            
            fprintf('\nThe following config fields were set to DEFAULTS and were NOT USED\n');
            sel = find(~used & ~original);
            if numel(sel)
              fprintf('  cfg.%s\n', key{sel});
            else
              fprintf('  <none>\n');
            end
          end % report
        end % report/cleanup
        
        if strcmp(cfg.trackconfig, 'cleanup')
          % remove the unused options from the configuration
          unusedkey = key(~used);
          for i=1:length(unusedkey)
            cfg = rmfield(cfg, unusedkey{i});
          end
        end
        
        % convert the configuration back to a struct
        cfg = struct(cfg);
      end
    end % off
    
  catch
    disp(lasterr);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check the size of fields in the cfg, remove large fields
% the max allowed size should be specified in cfg.checksize (this can be
% set with ft_defaults)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(checksize, 'yes') && ~isinf(cfg.checksize)
  cfg = checksizefun(cfg, cfg.checksize);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cfg] = checksizefun(cfg, max_size)

% first check the total size of the cfg
s = whos('cfg');
if (s.bytes <= max_size)
  return;
end

ignorefields = {'checksize', 'trl', 'trlold', 'event', 'artifact', 'artfctdef', 'previous'}; % these fields should never be removed
norecursion  = {'event'}; % these fields should not be handled recursively

fieldsorig = fieldnames(cfg);
for i=1:numel(fieldsorig)
  for k=1:numel(cfg)  % process each structure in a struct-array
    
    if any(strcmp(fieldsorig{i}, ignorefields))
      % keep this field, regardless of its size
      continue
      
    elseif iscell(cfg(k).(fieldsorig{i}))
      % run recursively on each struct element that is contained in the cell-array
      for j=1:numel(cfg(k).(fieldsorig{i}))
        if isstruct(cfg(k).(fieldsorig{i}){j})
          cfg(k).(fieldsorig{i}){j} = checksizefun(cfg(k).(fieldsorig{i}){j}, max_size);
        end
      end
      
    elseif isstruct(cfg(k).(fieldsorig{i})) && ~any(strcmp(fieldsorig{i}, norecursion))
      % run recursively on a struct field
      cfg(k).(fieldsorig{i}) = checksizefun(cfg(k).(fieldsorig{i}), max_size);
      
    else
      % determine the size of the field and remove it if too large
      temp = cfg(k).(fieldsorig{i});
      s = whos('temp');
      if s.bytes>max_size
        cfg(k).(fieldsorig{i}) = 'empty - this was cleared by checkconfig';
      end
      clear temp
      
    end
  end % for numel(cfg)
end % for each of the fieldsorig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION converts a cell array of structure arrays into a structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newNeighbours = fixneighbours(neighbours)
newNeighbours = struct;
for i=1:numel(neighbours)
  if i==1, newNeighbours = neighbours{i};    end;
  newNeighbours(i) = neighbours{i};
end
