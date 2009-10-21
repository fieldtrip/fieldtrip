function [cfg] = checkconfig(cfg, varargin)

% CHECKCONFIG checks the input cfg of the main FieldTrip functions.
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
%   [cfg] = checkconfig(cfg, ...)
%
% The behaviour of checkconfig can be controlled by the following cfg options,
% which can be set as global fieldtrip defaults (see FIELDTRIPDEFS):
%   cfg.checkconfig = 'pedantic', 'loose' or 'silent' (control the feedback behaviour of checkconfig)
%   cfg.trackconfig = 'cleanup', 'report' or 'off'
%   cfg.checksize   = number in bytes, can be inf (set max size allowed for output cfg fields)
%
% Optional input arguments should be specified as key-value pairs and can include
%   renamed         = {'old',  'new'}        % list the old and new option
%   renamedval      = {'opt',  'old', 'new'} % list option and old and new value
%   required        = {'opt1', 'opt2', etc.} % list the required options
%   deprecated      = {'opt1', 'opt2', etc.} % list the deprecated options
%   unused          = {'opt1', 'opt2', etc.} % list the unused options, these will be removed and a warning is issued
%   forbidden       = {'opt1', 'opt2', etc.} % list the forbidden options, these result in an error
%   createsubcfg    = {'subname', etc.}      % list the names of the subcfg
%   dataset2files   = 'yes', 'no'            % converts dataset into headerfile and datafile
%   checksize       = 'yes', 'no'            % remove large fields from the cfg
%   trackconfig     = 'on', 'off'            % start/end config tracking
%
% See also CHECKDATA, FIELDTRIPDEFS

% Copyright (C) 2007-2008, Robert Oostenveld, Saskia Haegens
%
% $Log: checkconfig.m,v $
% Revision 1.18  2009/08/05 13:03:55  roboos
% changed selection of file based on 'gui'
%
% Revision 1.17  2009/07/15 12:10:07  jansch
% added subspace and keepsubspace for subcfg dics and lcmv
%
% Revision 1.16  2009/06/04 13:41:33  marvger
% renamed mvlap case
%
% Revision 1.15  2009/05/22 13:20:18  marvger
% added case for mvlap method
%
% Revision 1.14  2009/05/14 18:54:39  roboos
% added sam for createsubcfg
%
% Revision 1.13  2009/04/03 08:09:31  jansch
% added denoise for preproc and subspace
%
% Revision 1.12  2009/03/12 17:10:38  roboos
% fixed bug in cfg.dataformat/headerformat (dataset2files) which applied to vhdr as input filename with a *.seg containing the data
%
% Revision 1.11  2009/02/04 09:09:59  roboos
% fixed filename to headerfile/datafile cvonversion in case of ctf_old
%
% Revision 1.10  2009/01/21 11:27:02  marvger
% automatically set dataformat and headerformat if unspecified
%
% Revision 1.9  2009/01/20 15:55:54  sashae
% cfg.trkcfgcount keeps track of number of times that trackconfig has been turned on/off,
% to prevent that nested functions turn off configtracking (i.e. only at end of the main function
% report/cleanup should be given)
%
% Revision 1.8  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.7  2008/12/16 15:37:23  sashae
% if cfg.checkconfig='silent' do not display report for trackconfig
%
% Revision 1.6  2008/12/04 19:11:11  sashae
% added silent/loose/pedantic feedback for 'renamed' and 'renamedval'
%
% Revision 1.5  2008/12/04 11:47:03  jansch
% added fixedori when beamformer_dics
%
% Revision 1.4  2008/12/02 17:51:28  sashae
% added artfctdef to 'ignorefields'
%
% Revision 1.3  2008/11/21 13:16:10  jansch
% added fixedori to be passed on in the case of lcmv (thanks to Joachim).
%
% Revision 1.2  2008/11/21 10:14:30  sashae
% added list of fields that should be ignored by trackconfig and checksize
%
% Revision 1.1  2008/11/13 09:55:36  roboos
% moved from fieldtrip/private, fileio or from roboos/misc to new location at fieldtrip/public
%
% Revision 1.21  2008/11/12 11:20:52  sashae
% change in configtracking: now ignores cfg.checksize
%
% Revision 1.20  2008/11/11 18:29:45  sashae
% some changes in configtracking
%
% Revision 1.19  2008/11/11 15:24:16  sashae
% updated documentation, improved checksize
%
% Revision 1.18  2008/11/11 13:12:59  roboos
% added/improved size checking
%
% Revision 1.17  2008/11/11 10:40:59  sashae
% first implementation of checksize: checks for large fields in the cfg and removes them
%
% Revision 1.16  2008/11/10 12:16:12  roboos
% clarified the handling of the configuration tracking options
%
% Revision 1.15  2008/11/06 15:19:23  sashae
% several updates in configuration tracking
%
% Revision 1.14  2008/11/02 10:56:34  roboos
% explicit sharing of code for dataset2files with fileio read_header/data
%
% Revision 1.13  2008/10/28 19:20:18  sashae
% configtracking can be turned on/off with trackconfig.
% user can request report on unused options with cfg.trackconfig='report', or get cleaned cfg with cfg.trackconfig='cleanup'.
%
% Revision 1.12  2008/10/13 13:38:59  sashae
% change in dataset2files code: empty dataset/headerfile/datafile fields are removed
%
% Revision 1.11  2008/10/13 12:41:33  jansch
% added projectmom for lcmv (in contrast to pcc the projection is done within
% beamformer_lcmv, instead of in sourcedescriptives). Probably this should
% change back when a clean version of sourcedescriptives is implemented.
%
% Revision 1.10  2008/10/10 12:33:04  sashae
% incorporated dataset2files
%
% Revision 1.9  2008/10/10 10:50:34  sashae
% updated documentation
%
% Revision 1.8  2008/10/02 14:06:21  roboos
% get the fields from ft_default and add them to the cfg structure
% implemented cfg.checkconfig=silent/loose/pedantic (default is in ft_default, i.e. fieldtripdefs function)
%
% Revision 1.7  2008/10/02 12:35:14  roboos
% added option "unused", renamed tracking to configtracking
%
% Revision 1.6  2008/10/01 15:45:29  sashae
% incorporated createsubcfg
%
% Revision 1.5  2008/09/30 13:05:19  roboos
% adedd first version of configuration tracking for testing
%
% Revision 1.4  2008/09/23 12:05:33  sashae
% some small changes; checkconfig can now handle empty cfgs
%
% Revision 1.3  2008/09/18 10:01:57  sashae
% added 'renamedval' which checks/adjusts renamed values
%
% Revision 1.2  2008/09/18 08:33:48  sashae
% new version: checks required, renamed, deprecated and forbidden configuration options,
% adjusts where possible and gives warning/error messages. to be implemented in all main
% fieldtrip functions, comparable to checkdata
%
% Revision 1.1  2008/07/08 15:39:22  roboos
% initial version for Saskia to work on
%

if isempty(cfg)
  cfg = struct; % ensure that it is an empty struct, not empty double
end

global ft_default
if isempty(ft_default)
  ft_default = struct;
end
fieldsused = fieldnames(ft_default);
for i=1:length(fieldsused)
  fn = fieldsused{i};
  if ~isfield(cfg, fn),
    cfg.(fn) = ft_default.(fn);
  end
end

renamed         = keyval('renamed',         varargin);
renamedval      = keyval('renamedval',      varargin);
required        = keyval('required',        varargin);
deprecated      = keyval('deprecated',      varargin);
unused          = keyval('unused',          varargin);
forbidden       = keyval('forbidden',       varargin);
createsubcfg    = keyval('createsubcfg',    varargin);
dataset2files   = keyval('dataset2files',   varargin);
checksize       = keyval('checksize',       varargin); if isempty(checksize), checksize = 'off';  end
trackconfig     = keyval('trackconfig',     varargin);

if ~isempty(trackconfig) && strcmp(trackconfig, 'on')
  % infer from the user configuration whether tracking should be enabled
  if isfield(cfg, 'trackconfig') && (strcmp(cfg.trackconfig, 'report') || strcmp(cfg.trackconfig, 'cleanup'))
    trackconfig = 'on'; % turn on configtracking if user requests report/cleanup
  else
    trackconfig = []; % disable configtracking if user doesn't request report/cleanup
  end
end

% these should be cell arrays and not strings
if ischar(required),   required   = {required};   end
if ischar(deprecated), deprecated = {deprecated}; end
if ischar(unused),     unused     = {unused};     end
if ischar(forbidden),  forbidden  = {forbidden};  end

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
  fieldsused = fieldnames(cfg);
  if any(strcmp(renamed{1}, fieldsused))
    cfg = setfield(cfg, renamed{2}, (getfield(cfg, renamed{1})));
    cfg = rmfield(cfg, renamed{1});
    if silent
      % don't mention it
    elseif loose
      warning(sprintf('use cfg.%s instead of cfg.%s', renamed{2}, renamed{1}));
    elseif pedantic
      error(sprintf('use cfg.%s instead of cfg.%s', renamed{2}, renamed{1}));
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rename old to new value, give warning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(renamedval) && isfield(cfg, renamedval{1})
  if strcmpi(getfield(cfg, renamedval{1}), renamedval{2})
    cfg = setfield(cfg, renamedval{1}, renamedval{3});
    if silent
      % don't mention it
    elseif loose
      warning(sprintf('use cfg.%s=%s instead of cfg.%s=%s', renamedval{1}, renamedval{3}, renamedval{1}, renamedval{2}));
    elseif pedantic
      error(sprintf('use cfg.%s=%s instead of cfg.%s=%s', renamedval{1}, renamedval{3}, renamedval{1}, renamedval{2}));
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
    error(sprintf('The field cfg.%s is required\n', required{ia}));
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
      warning(sprintf('The option cfg.%s is deprecated, support is no longer guaranteed\n', deprecated{ismember(deprecated, fieldsused)}));
    elseif pedantic
      error(sprintf('The option cfg.%s is deprecated, support is no longer guaranteed\n', deprecated{ismember(deprecated, fieldsused)}));
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
      warning(sprintf('The field cfg.%s is unused, it will be removed from your configuration\n', unused{ismember(unused, fieldsused)}));
    elseif pedantic
      error(sprintf('The field cfg.%s is unused\n', unused{ismember(unused, fieldsused)}));
    end
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
      warning(sprintf('The field cfg.%s is forbidden, it will be removed from your configuration\n', forbidden{ismember(forbidden, fieldsused)}));
    elseif pedantic
      error(sprintf('The field cfg.%s is forbidden\n', forbidden{ismember(forbidden, fieldsused)}));
    end
  end
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
          'blc'
          'blcwindow'
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
          'keepfilter'
          'keepmom'
          'lambda'
          'normalize'
          'normalizeparam'
          'powmethod'
          'projectnoise'
          'reducerank'
          'keepcsd'
          'realfilter'
	  'subspace'
	  'keepsubspace'
          };

      case 'lcmv'
        fieldname = {
          'feedback'
          'fixedori'
          'keepfilter'
          'keepmom'
          'lambda'
          'normalize'
          'normalizeparam'
          'powmethod'
          'projectnoise'
          'projectmom'
          'reducerank'
          'keepcov'
	  'subspace'
	  'keepsubspace'
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
          };

      case {'mne', 'loreta', 'rv'}
        fieldname = {
          'feedback'
          };

      case 'music'
        fieldname = {
          'feedback'
          'numcomponent'
          };
        
      case 'sam'
        fieldname = {
          'meansphereorigin'
          'spinning'
          'feedback'
          'lambda'
          'normalize'
          'normalizeparam'
          'reducerank'
          };

      case 'mvl'
        fieldname = {};
        
      otherwise
        error('unexpected name of the subfunction');
        fieldname = {};

    end % switch subname

    for i=1:length(fieldname)
      if ~isfield(subcfg, fieldname{i}) && isfield(cfg, fieldname{i})
        subcfg = setfield(subcfg, fieldname{i}, getfield(cfg, fieldname{i}));  % set it in the subconfiguration
        cfg = rmfield(cfg, fieldname{i});                                      % remove it from the main configuration
      end
    end

    % copy the substructure back into the main configuration structure
    cfg = setfield(cfg, subname, subcfg);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dataset2files
%
% Converts cfg.dataset into cfg.headerfile and cfg.datafile if neccessary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(dataset2files) && strcmp(dataset2files, 'yes')

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
    if strcmp(cfg.dataset, 'gui');
      [f, p] = uigetfile('*.*', 'Select a file');
      if isequal(f, 0)
        error('User pressed cancel');
      else
        d = fullfile(p, f);
      end
      cfg.dataset = d;
    end

    % the following code is shared with fileio read_header/read_data
    % therefore the three local variables are used outside of the cfg
    filename   = cfg.dataset;
    datafile   = [];
    headerfile = [];
    switch filetype(filename)
      case '4d_pdf'
        datafile   = filename;
        headerfile = [datafile '.m4d'];
        sensorfile = [datafile '.xyz'];
      case {'4d_m4d', '4d_xyz'}
        datafile   = filename(1:(end-4)); % remove the extension
        headerfile = [datafile '.m4d'];
        sensorfile = [datafile '.xyz'];
      case '4d'
        [path, file, ext] = fileparts(filename);
        datafile   = fullfile(path, [file,ext]);
        headerfile = fullfile(path, [file,ext]);
        configfile = fullfile(path, 'config');
      case {'ctf_ds', 'ctf_old'}
        % convert CTF filename into filenames
        [path, file, ext] = fileparts(filename);
        if any(strcmp(ext, {'.res4' '.meg4', '.1_meg4' '.2_meg4' '.3_meg4' '.4_meg4' '.5_meg4' '.6_meg4' '.7_meg4' '.8_meg4' '.9_meg4'}))
          filename = path;
          [path, file, ext] = fileparts(filename);
        end
        if isempty(path) && isempty(file)
          % this means that the dataset was specified as the present working directory, i.e. only with '.'
          filename = pwd;
          [path, file, ext] = fileparts(filename);
        end
        headerfile = fullfile(filename, [file '.res4']);
        datafile   = fullfile(filename, [file '.meg4']);
        if length(path)>3 && strcmp(path(end-2:end), '.ds')
          filename = path; % this is the *.ds directory
        end
      case 'ctf_meg4'
        [path, file, ext] = fileparts(filename);
        if isempty(path)
          path = pwd;
        end
        headerfile = fullfile(path, [file '.res4']);
        datafile   = fullfile(path, [file '.meg4']);
        if length(path)>3 && strcmp(path(end-2:end), '.ds')
          filename = path; % this is the *.ds directory
        end
      case 'ctf_res4'
        [path, file, ext] = fileparts(filename);
        if isempty(path)
          path = pwd;
        end
        headerfile = fullfile(path, [file '.res4']);
        datafile   = fullfile(path, [file '.meg4']);
        if length(path)>3 && strcmp(path(end-2:end), '.ds')
          filename = path; % this is the *.ds directory
        end
      case 'brainvision_vhdr'
        [path, file, ext] = fileparts(filename);
        headerfile = fullfile(path, [file '.vhdr']);
        if exist(fullfile(path, [file '.eeg']))
          datafile   = fullfile(path, [file '.eeg']);
        elseif exist(fullfile(path, [file '.seg']))
          datafile   = fullfile(path, [file '.seg']);
        elseif exist(fullfile(path, [file '.dat']))
          datafile   = fullfile(path, [file '.dat']);
        end
      case 'brainvision_eeg'
        [path, file, ext] = fileparts(filename);
        headerfile = fullfile(path, [file '.vhdr']);
        datafile   = fullfile(path, [file '.eeg']);
      case 'brainvision_seg'
        [path, file, ext] = fileparts(filename);
        headerfile = fullfile(path, [file '.vhdr']);
        datafile   = fullfile(path, [file '.seg']);
      case 'brainvision_dat'
        [path, file, ext] = fileparts(filename);
        headerfile = fullfile(path, [file '.vhdr']);
        datafile   = fullfile(path, [file '.dat']);
      case 'fcdc_matbin'
        [path, file, ext] = fileparts(filename);
        headerfile = fullfile(path, [file '.mat']);
        datafile   = fullfile(path, [file '.bin']);
      otherwise
        % convert filename into filenames, assume that the header and data are the same
        datafile   = filename;
        headerfile = filename;
    end
    % end sharing with fileio read_header/read_data
    % put everything back into the cfg
    cfg.dataset    = filename;
    cfg.datafile   = datafile;
    cfg.headerfile = headerfile;

    % fill dataformat if unspecified
    if ~isfield(cfg,'dataformat') || isempty(cfg.dataformat)
      cfg.dataformat = filetype(datafile);
    end

    % fill dataformat if unspecified
    if ~isfield(cfg,'headerformat') || isempty(cfg.headerformat)
      cfg.headerformat = filetype(headerfile);
    end
    
  elseif ~isempty(cfg.datafile) && isempty(cfg.headerfile);
    % assume that the datafile also contains the header
    cfg.headerfile = cfg.datafile;
  elseif isempty(cfg.datafile) && ~isempty(cfg.headerfile);
    % assume that the headerfile also contains the data
    cfg.datafile = cfg.headerfile;
  end
  % remove empty fields (otherwise a subsequent check on required fields doesn't make any sense)
  if isempty(cfg.dataset),    cfg=rmfield(cfg, 'dataset');    end
  if isempty(cfg.headerfile), cfg=rmfield(cfg, 'headerfile'); end
  if isempty(cfg.datafile),   cfg=rmfield(cfg, 'datafile');   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% configtracking
%
% switch configuration tracking on/off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(trackconfig)
  try
    if strcmp(trackconfig, 'on') && isa(cfg, 'struct')
      % turn ON configuration tracking
      cfg = config(cfg);
      % remember that configtracking has been turned on
      cfg.trkcfgcount = 1;
    elseif strcmp(trackconfig, 'on') && isa(cfg, 'config')
      % remember how many times configtracking has been turned on
      cfg.trkcfgcount = cfg.trkcfgcount+1; % count the 'ONs'
    end

    if strcmp(trackconfig, 'off') && isa(cfg, 'config')
      % turn OFF configuration tracking, optionally give report and/or cleanup
      cfg.trkcfgcount=cfg.trkcfgcount-1; % count(down) the 'OFFs'

      if cfg.trkcfgcount==0 % only proceed when number of 'ONs' matches number of 'OFFs'
        cfg=rmfield(cfg, 'trkcfgcount');

        if strcmp(cfg.trackconfig, 'report') || strcmp(cfg.trackconfig, 'cleanup')
          % gather information about the tracked results
          r = access(cfg, 'reference');
          o = access(cfg, 'original');

          key = fieldnames(cfg);
          key = key(:)';

          ignorefields = {'checksize', 'trl', 'trlold', 'event', 'artifact', 'artfctdef', 'previous'}; % these fields should never be removed!
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
% set with fieldtripdefs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(checksize, 'yes') && ~isinf(cfg.checksize)
  cfg = checksizefun(cfg, cfg.checksize);
end

function [cfg] = checksizefun(cfg, max_size)

ignorefields = {'checksize', 'trl', 'trlold', 'event', 'artifact', 'artfctdef', 'previous'}; % these fields should never be removed!

fieldsorig = fieldnames(cfg);
for i=1:numel(fieldsorig)
  for k=1:numel(cfg)
    if ~isstruct(cfg(k).(fieldsorig{i})) && ~any(strcmp(fieldsorig{i}, ignorefields))
      % find large fields and remove them from the cfg, skip fields that should be ignored
      temp = cfg(k).(fieldsorig{i});
      s = whos('temp');
      if s.bytes>max_size
        cfg(k).(fieldsorig{i}) = 'empty - this was cleared by checkconfig';
      end
      %%% cfg(k).(fieldsorig{i})=s.bytes; % remember the size of each field for debugging purposes
    elseif isstruct(cfg(k).(fieldsorig{i}));
      % run recursively on subfields that are structs
      cfg(k).(fieldsorig{i}) = checksizefun(cfg(k).(fieldsorig{i}), max_size);
    elseif iscell(cfg(k).(fieldsorig{i})) && strcmp(fieldsorig{i}, 'previous')
      % run recursively on 'previous' fields that are cells
      for j=1:numel(cfg(k).(fieldsorig{i}))
        cfg(k).(fieldsorig{i}){j} = checksizefun(cfg(k).(fieldsorig{i}){j}, max_size);
      end
    end
  end
end
