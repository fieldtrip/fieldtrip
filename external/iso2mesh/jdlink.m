function varargout = jdlink(uripath, varargin)
%
%    data = jdlink(uripath)
%       or
%    [data, fname, cachepath] = jdlink(uripath, 'param1', value1, ...)
%
%    Download linked data files from URLs and store those in cached folders
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        uripath: a single string or a cell array of strings, containing
%                the http:// or https:// links pointing to the linked
%                resources
%        'param'/value pairs: (optional) additional options are supported,
%                including
%            showlink: [1]: print URL or cached file; 0 do not print.
%            showsize: [1]: print the total size of the linked files; 0 do not print.
%            regex: a regular expression that is used to filter the URL
%                 cell array; only those matching the pattern are being
%                 downloaded; this has no effect to a single URL input
%
%    output:
%        data: a cell array storing the parsed data of each linked file
%        fname: a cell array listing the path to each locally cached files
%        cachepath: a cell array listing the cache search path orders
%
%    examples:
%        data = loadjson('https://neurojson.io:7777/openneuro/ds000001');
%        anatfiles = jsonpath(data, '$..anat.._DataLink_');
%        data = jdlink(anatfiles, 'regex', 'sub-0[12].*\.nii');
%        jsonpath(data, '$..Dim')
%
%    license:
%        BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details
%
% -- this function is part of JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
%

opt = varargin2struct(varargin{:});
opt.showlink = jsonopt('showlink', 1, opt);
opt.showsize = jsonopt('showsize', 1, opt);

if (iscell(uripath))
    if (isfield(opt, 'regex'))
        haspattern = cellfun(@(x) isempty(regexp(x, opt.regex, 'once')), uripath);
        uripath(haspattern) = [];
    end
    if (isfield(opt, 'showsize'))
        totalsize = 0;
        nosize = 0;
        for i = 1:length(uripath)
            filesize = regexp(uripath{i}, '&size=(\d+)', 'tokens');
            if (~isempty(filesize) && ~isempty(filesize{1}))
                totalsize = totalsize + str2double(filesize{1});
            else
                nosize = nosize + 1;
            end
        end
        fprintf('total %d links, %.0f bytes, %d files with unknown size\n', length(uripath), totalsize, nosize);
    end
    alloutput = cell(1, nargout);
    for i = 1:length(uripath)
        [newdata, fname, cachepath] = downloadlink(uripath{i}, opt);
        if (nargout > 0)
            alloutput{1}{end + 1} = newdata;
            if (nargout > 1)
                alloutput{2}{end + 1} = fname;
                if (nargout > 2)
                    alloutput{3}{end + 1} = cachepath;
                end
            end
        end
    end
    if (length(uripath) == 1)
        alloutput = cellfun(@(x) x{1}, alloutput, 'UniformOutput', false);
    end
    varargout = alloutput;
elseif (ischar(uripath) || isa(uripath, 'string'))
    [varargout{1:nargout}] = downloadlink(uripath, opt);
end

%%
function [newdata, fname, cachepath] = downloadlink(uripath, opt)
newdata = [];
[cachepath, filename] = jsoncache(uripath);
if (iscell(cachepath) && ~isempty(cachepath))
    if (opt.showlink)
        fprintf(1, 'downloading from URL: %s\n', uripath);
    end

    fname = [cachepath{1} filesep filename];
    fpath = fileparts(fname);
    if (~exist(fpath, 'dir'))
        mkdir(fpath);
    end
    if (exist('websave'))
        websave(fname, uripath);
    else
        rawdata = urlread(uripath);
        fid = fopen(fname, 'wb');
        if (fid == 0)
            error('can not save URL to cache at path %s', fname);
        end
        fwrite(fid, uint8(rawdata));
        fclose(fid);
    end
    newdata = loadjd(fname, opt);
elseif (~iscell(cachepath) && exist(cachepath, 'file'))
    if (opt.showlink)
        fprintf(1, 'loading from cache: %s\n', cachepath);
    end
    fname = cachepath;
    newdata = loadjd(fname, opt);
end
