function result=moxunit_fieldtrip_util_read_test_mfile_properties(fn)
% Read properties from Fieldtrip test file
%
% result=moxunit_fieldtrip_util_read_mfile_properties(fn)
%
% Inputs:
%   fn                      Filename with Fieldtrip test
%
% Output:
%   result                  struct with zero or more of the possible
%                           fields:
%                           'mem'       Maximum memory, in bytes
%                           'walltime'  Maximum time, in seconds
%                           'dccnpath'  true if the test uses 'dccnpath'
%
%
%
% #   For MOxUnit_fieldtrip's copyright information and license terms,   #
% #   see the COPYING file distributed with MOxUnit_fieldtrip.           #

    lines=read_lines_from_file(fn);

    % different properties we are looking for.
    % Each should have a field .pattern with a regular expression that
    % matches a single line. A field .parser sets a function handle that
    % parses the remaining part of the line. Alternatively
    %
    props=struct();

    props.walltime.prefix='^%\s*WALLTIME\s*';
    props.walltime.parser=@moxunit_fieldtrip_util_parse_walltime;

    props.mem.prefix='^%\s*MEM\s*';
    props.mem.parser=@moxunit_fieldtrip_util_parse_mem;

    props.dependency.prefix='^%\s*TEST\s*';
    props.dependency.parser=@moxunit_fieldtrip_util_parse_dependency;

    props.dccnpath.parser=@moxunit_fieldtrip_util_parse_dccnpath;

    % see which of these properties is defined in the file
    keys=fieldnames(props);

    result=struct();
    for k=1:numel(keys)
        key=keys{k};
        prop=props.(key);

        if isfield(prop,'prefix')
            value=parse_with_prefix(lines,prop.prefix,prop.parser);
        else
            value=parse_from_all_lines(lines,prop.parser);
        end

        if ~isempty(value)
            result.(key)=value;
        end
    end

function value=parse_with_prefix(lines,prefix,parser)
    value=[];
    match=regexp(lines,prefix,'end','once');
    idxs=find(~cellfun(@isempty,match));

    for k=1:numel(idxs)
        idx=idxs(k);
        remainder=lines{idx}((match{idx}+1):end);

        value=parser(remainder);
        if ~isempty(value)
            break;
        end
    end

function value=parse_from_all_lines(lines,parser)
    value=[];
    match=cellfun(parser,lines,'UniformOutput',false);
    msk=~cellfun(@isempty,match);

    if any(msk)
        value=true;
        return
    end


function lines=read_lines_from_file(fn)
    fid=fopen(fn);
    cleaner=onCleanup(@()fclose(fid));

    all_chars=fread(fid,inf,'*char');
    lines=regexp(all_chars',sprintf('\n'),'split');
