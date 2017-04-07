function obj=filter(obj, varargin)
% filter test cases based on execution time, memory use, or file presence
%
% obj=filter(obj, key, value)
%
% Inputs:
%   obj                 MoxUnitFieldTripTestSuite instance.
%   'maxwalltime',w     skip test if 'walltime' is set and greater than w
%   'maxmem,m           skip test if 'mem' is set and greater than m
%   'loadfile',lf       skip test if it used dccnpath and lf is false
%   'exclude_if_prefix_equals_failed',e
%                       skip test if the filename starts with 'failed'
%
% Output:
%   obj                 MoxUnitFieldTripTestSuite instance with the test
%                       replaced by a function handle test case
%                       that raises a skipped test exception. The number of
%                       tests after applying this function is unchanged.
%
% #   For MOxUnit_fieldtrip's copyright information and license terms,   #
% #   see the COPYING file distributed with MOxUnit_fieldtrip.           #

    % get filters based on exclusion criteria
    filters=get_filter_cell(varargin{:});

    % apply the filter
    obj=apply_filter_to_test_suite(obj,filters);



function suite=apply_filter_to_test_suite(suite,filters)
    for k=1:countTestNodes(suite)
        test_node=getTestNode(suite,k);

        fixed_test_node=[];

        if isa(test_node,'MOxUnitFieldTripTestSuite')
            % recursive call
            fixed_test_node=apply_filter_to_test_suite(test_node,filters);

        elseif isa(test_node,'MOxUnitFieldTripTestCase')
            [do_filter_out,reason]=apply_filters(filters,test_node);
            if do_filter_out
                % replace existing test node
                name=getName(test_node);
                location=getLocation(test_node);
                handle=@()moxunit_throw_test_skipped_exception(reason);
                fixed_test_node=MOxUnitFunctionHandleTestCase(name,...
                                                            location,...
                                                            handle);
            end
        end

        if ~isempty(fixed_test_node)
            % update test node
            suite=setTestNode(suite, k, fixed_test_node);
        end
    end


function [do_filter_out,reason]=apply_filters(filters,test_node)
    for k=1:2:numel(filters)
        key=filters{k};
        filter_func=filters{k+1};

        do_filter_out=filter_func(test_node);

        if do_filter_out
            reason=sprintf('''%s'' excludes %s (%s)',...
                          key,getName(test_node),getLocation(test_node));
            return;
        end
    end

    reason='';




function filters=get_filter_cell(varargin)
    narg=numel(varargin);

    filters=cell(1,narg);

    for k=1:2:narg
        key=varargin{k};
        raw_value=varargin{k+1};

        switch key
            case 'maxwalltime'
                value=from_number_or_parser(raw_value,...
                                @moxunit_fieldtrip_util_parse_walltime);
                exclude=@(tc) get(tc,'walltime')>value;
            case 'maxmem'
                value=from_number_or_parser(raw_value,...
                                @moxunit_fieldtrip_util_parse_mem);
                exclude=@(tc) get(tc,'mem')>value;

            case 'loadfile'
                value=istrue(raw_value);
                exclude=@(tc) isequal(get(tc,'dccnpath'),true) && ~value;

            case 'exclude_if_prefix_equals_failed'
                value=istrue(raw_value);
                exclude=@(tc) value && filename_starts_with('failed',tc);

            otherwise
                error('unsupported key %s', key);
        end

        filters{k}=key;
        filters{k+1}=exclude;
    end

function tf=filename_starts_with(prefix,test_case)
    full_path=getLocation(test_case);
    [unused,filename]=fileparts(full_path);

    pat=['^' regexptranslate('escape',prefix)];
    tf=~isempty(regexp(filename,pat,'once'));






function value=from_number_or_parser(raw_value, parser)
    if isnumeric(raw_value)
        value=raw_value;
    elseif ischar(raw_value)
        value=parser(raw_value);
    else
        error('illegal type %s', class(raw_value));
    end






