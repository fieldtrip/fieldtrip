function data = loadbidstsv(tsvfile, delim)
%
%    data = loadbidstsv(tsvfile)
%       or
%    data = loadbidstsv(tsvfile, delim)
%
%    Loading a BIDS-formatted .tsv (tab-separated values) or .tsv.gz file as a
%    struct; numerical fields are converted to floating-point data records
%    when possible; the header of the file is parsed to define sub-field
%    names
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        tsvfile: the path to the .tsv file
%        delim: (optional) if not set, tab ('\t') is used as column delimiter
%
%    examples:
%        data = loadbidstsv('participants.tsv');
%
%    license:
%        BSD license, see LICENSE_BSD.txt files for details
%
% -- this function is part of JBIDS toolbox (https://neurojson.org/#software)
%

if (nargin < 2)
    delim = sprintf('\t');
end

data = struct;

if (~isempty(regexp(tsvfile, '\.[Gg][Zz]$', 'once')))
    finput = fopen(tsvfile, 'rb');
    tsvdata = fread(finput, inf, 'uint8=>uint8');
    fclose(finput);

    if (~exist('gzipdecode', 'file'))
        error('To process zipped files, you must install gzipdecode.m from the JSONLab toolbox: http://github.com/NeuroJSON/jsonlab');
    end
    fid = char(gzipdecode(tsvdata));
    clear tsvdata;
    [header, endpos] = regexp(fid, '([^\n\r]*)', 'once', 'tokens', 'end');
    if (~isempty(header))
        header = header{1};
        fid = fid((endpos + 1):end);
    end
else
    fid = fopen(tsvfile, 'rb');
    header = fgetl(fid);
    header = regexprep(header, '\s*$', '');
end

if (isempty(header))
    return
end

if (exist('strsplit'))
    cols = strsplit(header, delim);
else
    cols = regexp(header, '\t*([^\t]*)\t*', 'tokens');
    cols = cellfun(@(x) x{:}, cols, 'uniformoutput', 0);
end
cols = cellfun(@encodevarname, cols, 'uniformoutput', 0);
if (~isempty(cols))
    body = textscan(fid, [repmat('%s\t', 1, length(cols) - 1), '%s'], 'delimiter', '\t');
    if (length(body) ~= length(cols))
        error('invalid tsv');
    end
    for i = 1:length(body)
        try
            % bodynum = cellfun(@(x) sscanf(regexprep(x,'^n/a$','NaN'), '%f'), body{i}, 'uniformoutput', 0);
            bodynum = cellfun(@(x) sscanf(x, '%f\t'), body{i}, 'uniformoutput', 0);
            if (exist('isna'))
                len = cellfun(@(x) numel(x) * (~isna(sum(x))), bodynum);
            else
                len = cellfun(@numel, bodynum);
            end
            if (any(len))
                body{i}(len > 0) = bodynum(len > 0);
                if (all(len))
                    body{i} = cell2mat(body{i});
                end
            end
        catch ME
            warning(ME.message);
        end
        data.(cols{i}) = body{i}(:).';
    end
end

fclose(fid);
