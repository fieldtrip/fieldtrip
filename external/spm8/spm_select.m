function varargout = spm_select(varargin)
% File selector
% FORMAT [t,sts] = spm_select(n,typ,mesg,sel,wd,filt,frames)
%     n    - Number of files
%            A single value or a range.  e.g.
%            1       - Select one file
%            Inf     - Select any number of files
%            [1 Inf] - Select 1 to Inf files
%            [0 1]   - select 0 or 1 files
%            [10 12] - select from 10 to 12 files
%     typ  - file type
%           'any'   - all files
%           'image' - Image files (".img" and ".nii")
%                     Note that it gives the option to select
%                     individual volumes of the images.
%           'mesh'  - Surface mesh files (".gii" or ".mat")
%           'xml'   - XML files
%           'mat'   - Matlab .mat files or .txt files (assumed to contain
%                     ASCII representation of a 2D-numeric array)
%           'batch' - SPM batch files (.m, .mat and XML)
%           'dir'   - select a directory
%           Other strings act as a filter to regexp.  This means
%           that e.g. DCM*.mat files should have a typ of '^DCM.*\.mat$'
%      mesg - a prompt (default 'Select files...')
%      sel  - list of already selected files
%      wd   - Directory to start off in
%      filt - value for user-editable filter (default '.*')
%      frames - Image frame numbers to include (default '1')
%
%      t    - selected files
%      sts  - status (1 means OK, 0 means window quit)
%
% FORMAT [t,ind] = spm_select('Filter',files,typ,filt,frames)
% filter the list of files (cell or char array) in the same way as the
% GUI would do. There is an additional typ 'extimage' which will match
% images with frame specifications, too. Also, there is a typ 'extdir',
% which will match canonicalised directory names. The 'frames' argument
% is currently ignored, i.e. image files will not be filtered out if
% their frame numbers do not match.
% t returns the filtered list (cell or char array, depending on input),
% ind an index array, such that t = files{ind}, or t = files(ind,:).
%
% FORMAT cpath = spm_select('CPath',path,cwd)
% function to canonicalise paths: Prepends cwd to relative paths, processes
% '..' & '.' directories embedded in path.
% path     - string matrix containing path name
% cwd      - current working directory [default '.']
% cpath    - conditioned paths, in same format as input path argument
%
% FORMAT [files,dirs] = spm_select('List',direc,filt)
% Returns files matching the filter (filt) and directories within direc
% direc    - directory to search
% filt     - filter to select files with (see regexp) e.g. '^w.*\.img$'
% files    - files matching 'filt' in directory 'direc'
% dirs     - subdirectories of 'direc'
%
% FORMAT [files,dirs] = spm_select('ExtList',direc,filt,frames)
% As above, but for selecting frames of 4D NIfTI files
% frames   - vector of frames to select (defaults to 1, if not specified)
%
% FORMAT [files,dirs] = spm_select('FPList',direc,filt)
% FORMAT [files,dirs] = spm_select('ExtFPList',direc,filt,frames)
% As above, but returns files with full paths (i.e. prefixes direc to each)
% FORMAT [files,dirs] = spm_select('FPListRec',direc,filt)
% FORMAT [files,dirs] = spm_select('ExtFPListRec',direc,filt,frames)
% As above, but returns files with full paths (i.e. prefixes direc to
% each) and searches through sub directories recursively.
%
% FORMAT spm_select('prevdirs',dir)
% Add directory dir to list of previous directories.
% FORMAT dirs = spm_select('prevdirs')
% Retrieve list of previous directories.
%____________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

if ~isdeployed && ~exist('cfg_getfile','file')
    addpath(fullfile(spm('dir'),'matlabbatch'));
end;
% cfg_getfile expects and returns cellstr arguments for multi-line strings
if nargin > 0 && ischar(varargin{1}) && strcmpi(varargin{1},'filter') && ischar(varargin{2})
    varargin{2} = cellstr(varargin{2});
end
[t sts] = cfg_getfile(varargin{:});
% cfg_getfile returns cell arrays, convert to char arrays
if nargin > 0 && ischar(varargin{1})
    switch lower(varargin{1})
        case 'filter',
            if ischar(varargin{2})
                t = char(t);
            end
        case {'list','fplist','extlist','extfplist'},
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
