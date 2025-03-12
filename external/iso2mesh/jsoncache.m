function [cachepath, filename] = jsoncache(dbname, docname, filename, domain)
%
% cachepaths=jsoncache()
% [cachepath, filename]=jsoncache(hyperlink)
% [cachepath, tf]=jsoncache(filename)
% cachepath=jsoncache(dbname, docname, filename, domain)
%
% return the JSON cache folder where _DataLink_ hyperlinked data files are downloaded
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% input:
%    hyperlink: if a single input is provided, the function check if it is
%               a hyperlink starting with http://, https:// or ftp://, if
%               so, it trys to extract the database name, document name and
%               file name using NeuroJSON's standard link format as
%
%    https://neurojson.org/io/stat.cgi?dbname=..&docname=..&file=..&size=..
%
%               if the URL does not follow the above format, a SHA-256 hash
%               will be computed based on the full URL to produce filename;
%               dbname is set as the first 2 letters of the hash and
%               docname is set to the 3rd/4th letters of the hash; the
%               domain name is also extracted from the URL; if the URL
%               contains the file's suffix, it is appended to the filename.
%
%               if the string does not contain a link, or the link starts
%               with file://, it is treated as a local file path
%    dbname: the name of the NeuroJSON database (must exist)
%    docname: the name of the NeuroJSON dataset document (must exist)
%    filename: the name of the data file - may contain a relative folder
%    domain: optional, if not given, 'default' is used; otherwise, user can
%            specify customized domain name
%
% output:
%    cachepaths: if the linked file is found in any of the cache folders,
%            this returns the full path of the found file as a string;
%            otherwise, this stores a cell array listing the searched cache
%            folders in the search order
%    tf: if a file is found in the cache folder, this returns true;
%            otherwise, this contains the extracted file name.
%
%    the cached data files will be searched in the following order
%
%    [pwd '/.neurojson']             | on all OSes
%    /home/USERNAME/.neurojson       | on all OSes (per-user)
%    /home/USERNAME/.cache/neurojson | if on Linux (per-user)
%    /var/cache/neurojson            | if on Linux (system wide)
%    /home/USERNAME/Library/neurojson| if on MacOS (per-user)
%    /Library/neurojson              | if on MacOS (system wide)
%    C:\ProgramData\neurojson        | if on Windows (system wide)
%
%    if a global variable NEUROJSON_CACHE is set in 'base', it will be
%    used instead of the above search paths
%
%
% example:
%    [cachepath, filename] = jsoncache('https://neurojson.org/io/stat.cgi?action=get&db=openneuro&doc=ds000001&file=sub-01/anat/sub-01_inplaneT2.nii.gz&size=669578')
%    [cachepath, filename] = jsoncache('https://raw.githubusercontent.com/NeuroJSON/jsonlab/master/examples/example1.json')
%    [cachepath, filename] = jsoncache('https://neurojson.io:7777/adhd200/Brown')
%    [cachepath, filename] = jsoncache('https://neurojson.io:7777/openneuro/ds003805')
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

pathname = getenv('HOME');
cachepath = {[pwd filesep '.neurojson']};
if (strcmp(pathname, pwd) == 0)
    cachepath{end + 1} = [pathname filesep '.neurojson'];
end

if (ispc)
    cachepath{end + 1} = [getenv('PROGRAMDATA') filesep 'neurojson'];
elseif (ismac)
    cachepath{end + 1} = [pathname '/Library/neurojson'];
    cachepath{end + 1} = '/Library/neurojson';
else
    cachepath{end + 1} = [pathname '/.cache/neurojson'];
    cachepath{end + 1} = '/var/cache/neurojson';
end

if (nargin < 4)
    domain = 'default';
end

if (nargin == 1)
    link = dbname;
    if (~isempty(regexp(link, '^file://', 'once')) || isempty(regexp(link, '://', 'once')))
        filename = regexprep(link, '^file://', '');
        if (exist(filename, 'file'))
            cachepath = filename;
            filename = true;
            return
        end
    else
        if (~isempty(regexp(link, '^https*://neurojson.org/io/', 'once')))
            domain = 'io';
        else
            newdomain = regexprep(regexp(link, '^(https*|ftp)://[^\/?#:]+', 'match', 'once'), '^(https*|ftp)://', '');
            if (~isempty(newdomain))
                domain = newdomain;
            end
        end
        dbname = regexp(link, '(?<=db=)[^&]+', 'match', 'once');
        docname = regexp(link, '(?<=doc=)[^&]+', 'match', 'once');
        filename = regexp(link, '(?<=file=)[^&]+', 'match', 'once');
        if (isempty(filename) && strcmp(domain, 'neurojson.io'))
            ref = regexp(link, '^(https*|ftp)://neurojson.io(:\d+)*(?<dbname>/[^\/]+)(?<docname>/[^\/]+)(?<filename>/[^\/?]+)*', 'names', 'once');
            if (~isempty(ref))
                if (~isempty(ref.dbname))
                    dbname = ref.dbname(2:end);
                end
                if (~isempty(ref.docname))
                    docname = ref.docname(2:end);
                end
                if (~isempty(ref.filename))
                    filename = ref.filename(2:end);
                elseif (~isempty(dbname))
                    if (~isempty(docname))
                        filename = [docname '.json'];
                    else
                        filename = [dbname '.json'];
                    end
                end
            end
        end
        if (isempty(filename))
            filename = jdatahash(link);
            suffix = regexp(link, '\.\w{1,5}(?=([#&].*)*$)', 'match', 'once');
            filename = [filename suffix];
            if (isempty(dbname))
                dbname = filename(1:2);
            end
            if (isempty(docname))
                docname = filename(3:4);
            end
        end
    end
end

p = getvarfrom({'caller', 'base'}, 'NEUROJSON_CACHE');

if (nargin == 0 || nargin == 1 || nargin >= 3)
    if (~isempty(p))
        cachepath = [{p}, cachepath{:}];
    elseif (exist('dbname', 'var') && exist('docname', 'var'))
        cachepath = cellfun(@(x) [x filesep domain filesep dbname filesep docname], cachepath, 'UniformOutput', false);
    end
    if (exist('filename', 'var') && ~isempty(filename))
        for i = 1:length(cachepath)
            if (exist([cachepath{i} filesep filename], 'file'))
                cachepath = [cachepath{i} filesep filename];
                filename = true;
                return
            end
        end
    elseif (exist('link', 'var'))
        [pathname, fname, fext] = fileparts(link);
        filename = [fname fext];
    end
    if (~isempty(p))
        cachepath(2) = [];
    else
        cachepath(1) = [];
    end
    return
end
