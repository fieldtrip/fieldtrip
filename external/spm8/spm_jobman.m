function varargout = spm_jobman(varargin)
% Main interface for SPM Batch System
% This function provides a compatibility layer between SPM and matlabbatch.
% It translates spm_jobman callbacks into matlabbatch callbacks and allows
% to edit and run SPM5 style batch jobs.
%
% FORMAT job_id = spm_jobman
%        job_id = spm_jobman('interactive')
%        job_id = spm_jobman('interactive',job)
%        job_id = spm_jobman('interactive',job,node)
%        job_id = spm_jobman('interactive','',node)
% Runs the user interface in interactive mode. The job_id can be used to
% manipulate this job in cfg_util. Note that changes to the job in cfg_util
% will not show up in cfg_ui unless 'Update View' is called.
%
% FORMAT output_list = spm_jobman('serial')
%        output_list = spm_jobman('serial',job[,'',   input1,...inputN])
%        output_list = spm_jobman('serial',job ,node[,input1,...inputN])
%        output_list = spm_jobman('serial',''  ,node[,input1,...inputN])
% Runs the user interface in serial mode. I job is not empty, then node
% is silently ignored. Inputs can be a list of arguments. These are
% passed on to the open inputs of the specified job/node. Each input should
% be suitable to be assigned to item.val{1}. For cfg_repeat/cfg_choice
% items, input should be a cell list of indices input{1}...input{k} into
% item.value. See cfg_util('filljob',...) for details.
%
% FORMAT output_list = spm_jobman('run',job)
%        output_list = spm_jobman('run_nogui',job)
% Runs a job without X11 (as long as there is no graphics output from the
% job itself). The matlabbatch system does not need graphics output to run
% a job.
%
%     node - indicates which part of the configuration is to be used.
%            For example, it could be 'jobs.spatial.coreg'.
%
%     job  - can be the name of a jobfile (as a .m, .mat or a .xml), a
%            cellstr of filenames, a 'jobs'/'matlabbatch' variable or a
%            cell of 'jobs'/'matlabbatch' variables loaded from a jobfile.
%
% Output_list is a cell array and contains the output arguments from each
% module in the job. The format and contents of these outputs is defined in
% the configuration of each module (.prog and .vout callbacks).
%
% FORMAT spm_jobman('initcfg')
% Initialise cfg_util configuration and set path accordingly.
%
% FORMAT jobs = spm_jobman('spm5tospm8',jobs)
% Takes a cell list of SPM5 job structures and returns SPM8 compatible versions.
%
% FORMAT job = spm_jobman('spm5tospm8bulk',jobfiles)
% Takes a cell string with SPM5 job filenames and saves them in SPM8
% compatible format. The new job files will be MATLAB .m files. Their
% filenames will be derived from the input filenames. To make sure they are
% valid MATLAB script names they will be processed with
% genvarname(filename) and have a '_spm8' string appended to their
% filename.
%
% FORMAT spm_jobman('help',node)
%        spm_jobman('help',node,width)
% Creates a cell array containing help information.  This is justified
% to be 'width' characters wide. e.g.
%     h = spm_jobman('help','jobs.spatial.coreg.estimate');
%     for i=1:numel(h),fprintf('%s\n',h{i}); end;
%
% not implemented: FORMAT spm_jobman('defaults')
% Runs the interactive defaults editor.
%
% FORMAT [tag, job] = spm_jobman('harvest', job_id|cfg_item|cfg_struct)
% Take the job with id job_id in cfg_util and extract what is
% needed to save it as a batch job (for experts only). If the argument is a
% cfg_item or cfg_struct tree, it will be harvested outside cfg_util. 
% tag - tag of the root node of the current job/cfg_item tree
% job - harvested data from the current job/cfg_item tree
%
% FORMAT spm_jobman('pulldown')
% Creates a pulldown 'TASKS' menu in the Graphics window.
%
% not implemented: FORMAT spm_jobman('jobhelp')
% Creates a cell array containing help information specific for a certain
% job. Help is only printed for items where job specific help is
% present. This can be used together with spm_jobman('help') to create a
% job specific manual. This feature is available only on MATLAB R14SP2
% and higher.
%
% not implemented: FORMAT spm_jobman('chmod')
% Changes the modality for the TASKS pulldown.
%
% This code is based on earlier versions by John Ashburner, Philippe 
% Ciuciu and Guillaume Flandin.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% Copyright (C) 2008 Freiburg Brain Imaging

% Volkmar Glauche
% $Id$


if nargin==0
    h = cfg_ui;
    if nargout > 0, varargout = {h}; end
else
    cmd = lower(varargin{1});
    if any(strcmp(cmd, {'serial','interactive','run','run_nogui'}))
        if nargin > 1
            % sort out job/node arguments for interactive, serial, run cmds
            if nargin>=2 && ~isempty(varargin{2})
                % do not consider node if job is given
                if ischar(varargin{2}) || iscellstr(varargin{2})
                    jobs = load_jobs(varargin{2});
                elseif iscell(varargin{2})
                    if iscell(varargin{2}{1})
                        % assume varargin{2} is a cell of jobs
                        jobs = varargin{2};
                    else
                        % assume varargin{2} is a single job
                        jobs{1} = varargin{2};
                    end;
                end;
                [mljob comp] = canonicalise_job(jobs);
            elseif any(strcmp(cmd, {'interactive','serial'})) && nargin>=3 && isempty(varargin{2})
                % Node spec only allowed for 'interactive', 'serial'
                arg3       = regexprep(varargin{3},'^spmjobs\.','spm.');
                mod_cfg_id = cfg_util('tag2mod_cfg_id',arg3);
            else
                error('spm:spm_jobman:WrongUI', ...
                    'Don''t know how to handle this ''%s'' call.', lower(varargin{1}));
            end;
        end
    end;
    switch cmd
        case 'help'
            if (nargin < 2) || isempty(varargin{2})
                node = 'spm';
            else
                node = regexprep(varargin{2},'^spmjobs\.','spm.');
            end;
            if nargin < 3
                width = 60;
            else
                width = varargin{3};
            end;
            varargout{1} = cfg_util('showdocwidth', width, node);
            
        case 'initcfg'
            if ~isdeployed
                addpath(fullfile(spm('Dir'),'matlabbatch'));
                addpath(fullfile(spm('Dir'),'config'));
            end
            cfg_get_defaults('cfg_util.genscript_run', @genscript_run);
            cfg_util('initcfg'); % This must be the first call to cfg_util
            if ~spm('cmdline')
                f = cfg_ui('Visible','off'); % Create invisible batch ui
                f0 = findobj(f, 'Tag','MenuFile'); % Add entries to file menu
                f2 = uimenu(f0,'Label','Load SPM5 job', 'Callback',@load_job, ...
                    'HandleVisibility','off', 'tag','jobs', ...
                    'Separator','on');
                f3 = uimenu(f0,'Label','Bulk Convert SPM5 job(s)', ...
                    'Callback',@conv_jobs, ...
                    'HandleVisibility','off', 'tag','jobs');
            end
        case 'interactive',
            if exist('mljob', 'var')
                cjob = cfg_util('initjob', mljob);
            elseif exist('mod_cfg_id', 'var')
                if isempty(mod_cfg_id)
                    arg3 = regexprep(varargin{3},'^spmjobs\.','spm.');
                    warning('spm:spm_jobman:NodeNotFound', ...
                        ['Can not find executable node ''%s'' - running '...
                        'matlabbatch without default node.'], arg3);
                    cjob = cfg_util('initjob');
                else
                    cjob = cfg_util('initjob');
                    mod_job_id = cfg_util('addtojob', cjob, mod_cfg_id);
                    cfg_util('harvest', cjob, mod_job_id);
                end;
            else
                cjob = cfg_util('initjob');
            end;
            cfg_ui('local_showjob', findobj(0,'tag','cfg_ui'), cjob);
            if nargout > 0
                varargout{1} = cjob;
            end;
            
        case 'serial',
            if exist('mljob', 'var')
                cjob = cfg_util('initjob', mljob);
            else
                cjob = cfg_util('initjob');
                if nargin > 2
                    arg3 = regexprep(varargin{3},'^spmjobs\.','spm.');
                    [mod_cfg_id, item_mod_id] = cfg_util('tag2cfg_id', lower(arg3));
                    cfg_util('addtojob', cjob, mod_cfg_id);
                end;
            end;
            sts  = cfg_util('filljobui', cjob, @serial_ui, varargin{4:end});
            if sts
                cfg_util('run', cjob);
                if nargout > 0
                    varargout{1} = cfg_util('getalloutputs', cjob);
                end
            end;
            cfg_util('deljob', cjob);

        case {'run','run_nogui'}
            cjob = cfg_util('initjob', mljob);
            cfg_util('run', cjob);
            if nargout > 0
                varargout{1} = cfg_util('getalloutputs', cjob);
            end
            cfg_util('deljob', cjob);
            
        case {'spm5tospm8'}
            varargout{1} = canonicalise_job(varargin{2});

        case {'spm5tospm8bulk'}
            conv_jobs(varargin{2});
            
        case {'defaults'},
            warning('spm:spm_jobman:NotImplemented', 'Not yet implemented.');
            
        case {'pulldown'}
            pulldown;

        case {'chmod'}
            warning('spm:spm_jobman:NotImplemented', 'Callback ''%s'' not implemented.', varargin{1});

        case {'help'}
            warning('spm:spm_jobman:NotImplemented', 'Not yet implemented.');

        case {'jobhelp'}
            warning('spm:spm_jobman:NotImplemented', 'Callback ''%s'' not implemented.', varargin{1});

        case {'harvest'}
            if nargin == 1
                error('spm:spm_jobman:CantHarvest', ...
                        ['Can not harvest job without job_id. Please use ' ...
                         'spm_jobman(''harvest'', job_id).']);
            elseif cfg_util('isjob_id', varargin{2})
                [tag job] = cfg_util('harvest', varargin{2});
            elseif isa(varargin{2}, 'cfg_item')
                [tag job] = harvest(varargin{2}, varargin{2}, false, false);
            elseif isstruct(varargin{2})
                % try to convert into class before harvesting
                c = cfg_struct2cfg(varargin{2});
                [tag job] = harvest(c,c,false,false);
            else
                error('spm:spm_jobman:CantHarvestThis', ['Can not harvest ' ...
                                    'this argument.']);
            end;
            varargout{1} = tag;
            varargout{2} = job;

        otherwise
            error(['"' varargin{1} '" - unknown option']);
    end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [mljob, comp] = canonicalise_job(job)
% job: a cell list of job data structures.
% Check whether job is a SPM5 or matlabbatch job. In the first case, all
% items in job{:} should have a fieldname of either 'temporal', 'spatial',
% 'stats', 'tools' or 'util'. If this is the case, then job will be
% assigned to mljob{1}.spm, which is the tag of the SPM root
% configuration item.

comp = true(size(job));
mljob = cell(size(job));
for cj = 1:numel(job)
    for k = 1:numel(job{cj})
        comp(cj) = comp(cj) && any(strcmp(fieldnames(job{cj}{k}), ...
            {'temporal', 'spatial', 'stats', 'tools', 'util'}));
        if ~comp(cj)
            break;
        end;
    end;
    if comp(cj)
        tmp = convert_jobs(job{cj});
        for i=1:numel(tmp),
            mljob{cj}{i}.spm = tmp{i};
        end
    else
        mljob{cj} = job{cj};
    end;
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function conv_jobs(varargin)
% Select a list of jobs, canonicalise each of it and save as a .m file
% using gencode.
spm('pointer','watch');
if nargin == 0 || ~iscellstr(varargin{1})
    [fname sts] = spm_select([1 Inf], 'batch', 'Select job file(s)');
    fname = cellstr(fname);
    if ~sts, return; end;
else
    fname = varargin{1};
end;

joblist = load_jobs(fname);
for k = 1:numel(fname)
    if ~isempty(joblist{k})
        [p n e v] = spm_fileparts(fname{k});
        % Save new job as genvarname(*_spm8).m
        newfname = fullfile(p, sprintf('%s.m', ...
            genvarname(sprintf('%s_spm8', n))));
        fprintf('SPM5 job: %s\nSPM8 job: %s\n', fname{k}, newfname);
        cjob = cfg_util('initjob', canonicalise_job(joblist(k)));
        cfg_util('savejob', cjob, newfname);
        cfg_util('deljob', cjob);
    end;
end;
spm('pointer','arrow');
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function load_job(varargin)
% Select a single job file, canonicalise it and display it in GUI
[fname sts] = spm_select([1 Inf], 'batch', 'Select job file');
if ~sts, return; end;

spm('pointer','watch');
joblist = load_jobs(fname);
if ~isempty(joblist{1})
    spm_jobman('interactive',joblist{1});
end;
spm('pointer','arrow');
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function newjobs = load_jobs(job)
% Load a list of possible job files, return a cell list of jobs. Jobs can
% be either SPM5 (i.e. containing a 'jobs' variable) or SPM8/matlabbatch
% jobs. If a job file failed to load, an empty cell is returned in the
% list.
if ischar(job)
    filenames = cellstr(job);
else
    filenames = job;
end;
newjobs = {};
for cf = 1:numel(filenames)
    [p,nam,ext] = fileparts(filenames{cf});
    switch ext
        case '.xml',
            spm('Pointer','Watch');
            try
                loadxml(filenames{cf},'jobs');
            catch
                try
                    loadxml(filenames{cf},'matlabbatch');
                catch
                    warning('spm:spm_jobman:LoadFailed','LoadXML failed: ''%s''',filenames{cf});
                end;
            end;
            spm('Pointer');
        case '.mat'
            try
                S=load(filenames{cf});
                if isfield(S,'matlabbatch')
                    matlabbatch = S.matlabbatch;
                elseif isfield(S,'jobs')
                    jobs = S.jobs;
                else
                    warning('spm:spm_jobman:JobNotFound','No SPM5/SPM8 job found in ''%s''', filenames{cf});
                end
            catch
                warning('spm:spm_jobman:LoadFailed','Load failed: ''%s''',filenames{cf});
            end;
        case '.m'
            try
                fid = fopen(filenames{cf},'rt');
                str = fread(fid,'*char');
                fclose(fid);
                eval(str);
            catch
                warning('spm:spm_jobman:LoadFailed','Load failed: ''%s''',filenames{cf});
            end;
            if ~(exist('jobs','var') || exist('matlabbatch','var'))
                warning('spm:spm_jobman:JobNotFound','No SPM5/SPM8 job found in ''%s''', filenames{cf});
            end;
        otherwise
            warning('Unknown extension: ''%s''', filenames{cf});
    end;
    if exist('jobs','var')
        newjobs = [newjobs(:); {jobs}];
        clear jobs;
    elseif exist('matlabbatch','var')
        newjobs = [newjobs(:); {matlabbatch}];
        clear matlabbatch;
    end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function njobs = convert_jobs(jobs)
decel    = struct('spatial',struct('realign',[],'coreg',[],'normalise',[]),...
                 'temporal',[],...
                 'stats',[],...
                 'meeg',[],...
                 'util',[],...
                 'tools',struct('dartel',[]));
njobs  = {};
for i0 = 1:numel(jobs),
    tmp0  = fieldnames(jobs{i0});
    tmp0  = tmp0{1};
    if any(strcmp(tmp0,fieldnames(decel))),
        for i1=1:numel(jobs{i0}.(tmp0)),
            tmp1  = fieldnames(jobs{i0}.(tmp0){i1});
            tmp1  = tmp1{1};
            if ~isempty(decel.(tmp0)),
                if any(strcmp(tmp1,fieldnames(decel.(tmp0)))),
                    for i2=1:numel(jobs{i0}.(tmp0){i1}.(tmp1)),
                        njobs{end+1} = struct(tmp0,struct(tmp1,jobs{i0}.(tmp0){i1}.(tmp1){i2}));
                    end
                else
                    njobs{end+1} = struct(tmp0,jobs{i0}.(tmp0){i1});
                end
            else
                njobs{end+1} = struct(tmp0,jobs{i0}.(tmp0){i1});
            end
        end
    else
        njobs{end+1} = jobs{i0};
    end
end
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function pulldown
fg = spm_figure('findwin','Graphics');
if isempty(fg), return; end;
set(0,'ShowHiddenHandles','on');
delete(findobj(fg,'tag','jobs'));
set(0,'ShowHiddenHandles','off');
f0 = uimenu(fg,'Label','TASKS', ...
            'HandleVisibility','off', 'tag','jobs');
f1 = uimenu(f0,'Label','BATCH', 'Callback',@cfg_ui, ...
            'HandleVisibility','off', 'tag','jobs');
f4 = uimenu(f0,'Label','SPM (interactive)', ...
            'HandleVisibility','off', 'tag','jobs', 'Separator','on');
cfg_ui('local_setmenu', f4, cfg_util('tag2cfg_id', 'spm'), ...
       @local_init_interactive, false);
f5 = uimenu(f0,'Label','SPM (serial)', ...
            'HandleVisibility','off', 'tag','jobs');
cfg_ui('local_setmenu', f5, cfg_util('tag2cfg_id', 'spm'), ...
       @local_init_serial, false);
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function local_init_interactive(varargin)
cjob = cfg_util('initjob');
mod_cfg_id = get(gcbo,'userdata');
cfg_util('addtojob', cjob, mod_cfg_id);
cfg_ui('local_showjob', findobj(0,'tag','cfg_ui'), cjob);
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function local_init_serial(varargin)
mod_cfg_id = get(gcbo,'userdata');
cjob = cfg_util('initjob');
cfg_util('addtojob', cjob, mod_cfg_id);
sts = cfg_util('filljobui', cjob, @serial_ui);
if sts
    cfg_util('run', cjob);
end;
cfg_util('deljob', cjob);
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [val sts] = serial_ui(item)
% wrapper function to translate cfg_util('filljobui'... input requests into
% spm_input/cfg_select calls.
sts = true;
switch class(item),
    case 'cfg_choice',
        labels = cell(size(item.values));
        values = cell(size(item.values));
        for k = 1:numel(item.values)
            labels{k} = item.values{k}.name;
            values{k} = k;
        end;
        val = spm_input(item.name, 1, 'm', labels, values);
    case 'cfg_menu',
        val = spm_input(item.name, 1, 'm', item.labels, item.values);
        val = val{1};
    case 'cfg_repeat',
        labels = cell(size(item.values));
        values = cell(size(item.values));
        for k = 1:numel(item.values)
            labels{k} = item.values{k}.name;
            values{k} = k;
        end;
        % enter at least item.num(1) values
        for k = 1:item.num(1)
            val(k) = spm_input(sprintf('%s(%d)', item.name, k), 1, 'm', ...
                               labels, values);
        end;
        % enter more (up to varargin{3}(2) values
        labels = {labels{:} 'Done'};
        % values is a cell list of natural numbers, use -1 for Done
        values = {values{:} -1}; 
        while numel(val) < item.num(2)
            val1 = spm_input(sprintf('%s(%d)', item.name, numel(val)+1), 1, ...
                             'm', labels, values);
            if val1{1} == -1
                break;
            else
                val(end+1) = val1;
            end;
        end;
    case 'cfg_entry',
        val = spm_input(item.name, 1, item.strtype, '', item.num, ...
                        item.extras);
    case 'cfg_files',
        [t,sts] = cfg_getfile(item.num, item.filter, item.name, '', ...
                              item.dir, item.ufilter);
        if sts
            val = cellstr(t);
        else
            val = {};
            error('File selector was closed.');
        end;
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [code cont] = genscript_run
% Return code snippet to initialise SPM defaults and run a job generated by
% cfg_util('genscript',...) through spm_jobman.
modality = spm('CheckModality');
code{1}  = sprintf('spm(''defaults'', ''%s'');', modality);
code{2}  = 'spm_jobman(''serial'', jobs, '''', inputs{:});';
cont     = false;
