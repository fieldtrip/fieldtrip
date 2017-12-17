function varargout = spm_select(varargin)
% File selector
% FORMAT [t,sts] = spm_select(n,typ,mesg,sel,wd,filt,frames)
% n      - number of files [Default: Inf]
%          A single value or a range.  e.g.
%          1       - select one file
%          Inf     - select any number of files
%          [1 Inf] - select 1 to Inf files
%          [0 1]   - select 0 or 1 files
%          [10 12] - select from 10 to 12 files
% typ    - file type [Default: 'any']
%          'any'   - all files
%          'image' - Image files (".img" and ".nii")
%                    Note that it gives the option to select individuals
%                    volumes of the images.
%          'mesh'  - Surface mesh files (".gii")
%          'xml'   - XML files
%          'mat'   - MATLAB .mat files or .txt files (assumed to contain
%                    ASCII representation of a 2D-numeric array)
%          'batch' - SPM batch files (.m or .mat)
%          'dir'   - select a directory
%          Other strings act as a filter to regexp. This means that
%          e.g. DCM*.mat files should have a typ of '^DCM.*\.mat$'
% mesg   - a prompt [Default: 'Select files...']
% sel    - list of already selected files [Default: {}]
% wd     - directory to start off in [Default: pwd]
% filt   - value for user-editable filter [Default: '.*']
% frames - image frame numbers to include [Default: '1']
%
% t      - selected files
% sts    - status (1 means OK, 0 means window quit)
%
% FORMAT [files,dirs] = spm_select('List',direc,filt)
% Return files matching the filter 'filt' and directories within 'direc'
% direc  - directory to search [Default: pwd]
% filt   - filter to select files with regexp, e.g. '^w.*\.img$' [Default: '.*']
%
% files  - files matching 'filt' in directory 'direc'
% dirs   - subdirectories of 'direc'
%
% FORMAT [files,dirs] = spm_select('ExtList',direc,filt,frames)
% As above, but for selecting frames of 4D NIfTI files
% frames - vector of frames to select (defaults to 1, if not specified).
%          If the frame number is Inf, all frames for the matching images
%          are listed. 
%
% FORMAT [files,dirs] = spm_select('FPList',direc,filt)
% FORMAT [files,dirs] = spm_select('ExtFPList',direc,filt,frames)
% As above, but return files with full paths (i.e. prefixes 'direc' to each)
%
% FORMAT [files,dirs] = spm_select('FPListRec',direc,filt)
% FORMAT [files,dirs] = spm_select('ExtFPListRec',direc,filt,frames)
% As above, but return files with full paths (i.e. prefixes 'direc' to each)
% and search through sub directories recursively.
%
% FORMAT [dirs] = spm_select('List',direc,'dir',filt)
% FORMAT [dirs] = spm_select('FPList',direc,'dir',filt)
% FORMAT [dirs] = spm_select('FPListRec',direc,'dir',filt)
% Return directory names matching filter 'filt' within 'direc'
%__________________________________________________________________________
% Copyright (C) 2005-2015 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_select.m 6530 2015-08-21 14:43:52Z guillaume $


% For developers:
%--------------------------------------------------------------------------
% FORMAT cpath = spm_select('CPath',path,cwd)
% Canonicalise paths: prepend cwd to relative paths, process '..' & '.'
% directories embedded in path.
% path   - string matrix containing path name
% cwd    - current working directory [Default: '.']
% cpath  - canonicalised paths, in same format as input path argument
%
% FORMAT [t,ind] = spm_select('Filter',files,typ,filt,frames)
% Filter the list of files (cell or char array) in the same way as the
% GUI would do. There is an additional typ 'extimage' which will match
% images with frame specifications, too. Also, there is a typ 'extdir',
% which will match canonicalised directory names. The 'frames' argument
% is currently ignored, i.e. image files will not be filtered out if
% their frame numbers do not match.
% t returns the filtered list (cell or char array, depending on input),
% ind an index array, such that t = files{ind}, or t = files(ind,:).
%
% FORMAT spm_select('PrevDirs',dir)
% Add directory dir to list of previous directories.
% FORMAT dirs = spm_select('PrevDirs')
% Retrieve list of previous directories.
%
% FORMAT files = spm_select('Expand',files)
% Return a list of image filenames with appended frame numbers.

persistent isInitSelect;
if isempty(isInitSelect)
    isInitSelect = true;
    spm_select('init');
    if nargin == 1 && strcmpi(varargin{1},'init'), return; end
end

%-Commands that are not passed to cfg_getfile
local_cmds = {'regfilter', 'init', 'expand'};

if nargin && ischar(varargin{1}) && any(strcmpi(varargin{1},local_cmds))
    switch lower(varargin{1})
        case 'init'
            if ~isdeployed && ~exist('cfg_getfile','file')
                addpath(fullfile(spm('dir'),'matlabbatch'));
            end
            spm_select('regfilter');
            spm_select('prevdirs',spm('Dir'));
            
        case 'regfilter'
            % Regexp based filters without special handlers
            cfg_getfile('regfilter', 'mesh',  {'\.gii$'});
            cfg_getfile('regfilter', 'gifti', {'\.gii$'});
            cfg_getfile('regfilter', 'nifti', {'\.nii$','\.img$'});
            % Filter for 3D images that handles frame expansion
            frames = cfg_entry;
            frames.name    = 'Frames';
            frames.tag     = 'frames';
            frames.strtype = 'n';
            frames.num     = [1 Inf];
            frames.val     = {1};
            cfg_getfile('regfilter', 'image',...
                {'.*\.nii(,\d+){0,2}$','.*\.img(,\d+){0,2}$'},...
                false, @spm_select_image, {frames});
            
        case 'expand'
            varargout{1} = spm_select_expand(varargin{2});
            if ischar(varargin{2}), varargout{1} = char(varargout{1}); end
    end
else
    needchar = false;
    % cfg_getfile expects cellstr arguments for multi-line strings
    if nargin > 1 && ischar(varargin{1}) && ...
            ismember(lower(varargin{1}),{'filter','cpath'}) && ischar(varargin{2})
        varargin{2} = cellstr(varargin{2});
        needchar = true;
    elseif nargin > 0 && ischar(varargin{1}) && ...
            ismember(lower(varargin{1}),{'extlist','extfplist','extfplistrec'})
        varargin{1} = varargin{1}(4:end);
        if nargin > 3
            varargin{5} = struct('frames', varargin{4});
        else
            varargin{5} = struct('frames', Inf);
        end
        varargin{4} = varargin{3};
        varargin{3} = 'image';
    elseif nargin > 6 && isnumeric(varargin{1})
        varargin{7} = struct('frames', varargin{7});
    end
    
    [t, sts] = cfg_getfile(varargin{:});
    
    % cfg_getfile returns cellstr arrays, convert to char arrays
    if nargin > 0 && ischar(varargin{1})
        switch lower(varargin{1})
            case {'filter','cpath'}
                if needchar
                    t = char(t);
                end
            case {'list','fplist','extlist','extfplist','fplistrec','extfplistrec'}
                t = char(t);
                sts = char(sts);
        end
    else
        t = char(t);
    end

    varargout{1} = t;
    if nargout > 1
        varargout{2} = sts;
    end
end


%==========================================================================
% FUNCTION varargout = spm_select_image(cmd, varargin)
%==========================================================================
function varargout = spm_select_image(cmd, varargin)
% Implements extended filtering for NIfTI images (including 3D frame
% selection)
switch lower(cmd)
    case 'list'
        dr     = varargin{1};
        files  = varargin{2};
        prms   = varargin{3};
        frames = prms.frames;
        ii = cell(1,numel(files));
        if numel(frames)==1 && isnan(frames)
            [ii{:}] = deal(1);
        elseif numel(frames)~=1 || frames(1)~=1
            % if domsg
            %     msg(ob,['Reading headers of ' num2str(numel(f)) ' images...']);
            % end
            for i=1:numel(files)
                try
                    n = spm_select_get_nbframes(fullfile(dr,files{i}));
                    d4 = (1:n)';
                catch
                    d4 = 1;
                end
                if all(isfinite(frames))
                    ii{i} = intersect(d4, frames(:))';
                else
                    ii{i} = d4(:)';
                end
            end
        elseif numel(frames)==1 && frames(1)==1
            [ii{:}] = deal(1);
        end

        % if domsg
        %     msg(ob,['Listing ' num2str(numel(f)) ' files...']);
        % end

        % Combine filename and frame number(s)
        nii      = cellfun(@numel, ii);
        cfiles   = cell(sum(nii),1);
        fi       = cell(numel(files),1);
        for k = 1:numel(fi)
            fi{k} = k*ones(1,nii(k));
        end
        ii = [ii{:}];
        fi = [fi{:}]';
        if numel(frames)==1 && isnan(frames)
            cfiles = files;
        else
            for i=1:numel(cfiles)
                cfiles{i} = sprintf('%s,%d', files{fi(i)}, ii(i));
            end
        end
        varargout{1} = cfiles;
        varargout{2} = fi;
    case 'filter'
        % Do not filter for frame numbers
        varargout{1} = varargin{1};
        varargout{2} = 1:numel(varargout{1});
end


%==========================================================================
% FUNCTION ofiles = spm_select_expand(ifiles)
%==========================================================================
function ofiles = spm_select_expand(ifiles)
ifiles = cellstr(ifiles);
ofiles = {};
for i=1:numel(ifiles)
    [p,f,e,n] = spm_fileparts(ifiles{i});
    if ~isempty(n)
        ofiles = [ofiles; ifiles{i}];
    else
        n = spm_select_get_nbframes(ifiles{i});
        vfiles = cell(n,1);
        for j=1:n
            vfiles{j} = [ifiles{i} ',' num2str(j)];
        end
        ofiles = [ofiles; vfiles];
    end
end


%==========================================================================
% FUNCTION n = spm_select_get_nbframes(file)
%==========================================================================
function n = spm_select_get_nbframes(file)
N   = nifti(file);
dim = [N.dat.dim 1 1 1 1 1];
n   = dim(4);

% % A faster, direct implementation
% [p,f,e] = fileparts(file);
% unchanged = true;
% ind = find(e==',');
% if ~isempty(ind)
%     unchanged = false;
%     e = e(1:(ind(1)-1));
% end
% switch e
%     case {'.img'}
%         unchanged = false;
%         e = '.hdr';
%     case {'.IMG'}
%         unchanged = false;
%         e = '.HDR';
% end
% if ~unchanged
%     file = fullfile(p,[f e]);
% end
% 
% n = [];
% 
% fp  = fopen(file,'r','native');
% if fp==-1, return; end
% 
% %sts = fseek(fp,344,'bof');
% %if sts==-1, fclose(fp); return; end
% %mgc = deblank(char(fread(fp,4,'uint8')'));
% 
% sts = fseek(fp,40,'bof');
% if sts==-1, fclose(fp); return; end
% dim = fread(fp,8,'*int16');
% fclose(fp);
% if isempty(dim), return; end
% if dim(1)<1 || dim(1)>7
%     dim = swapbytes(dim);
% end
% dim = double(dim(2:(dim(1)+1)))';
% dim = [dim 1 1 1 1 1];
% n   = dim(4);
