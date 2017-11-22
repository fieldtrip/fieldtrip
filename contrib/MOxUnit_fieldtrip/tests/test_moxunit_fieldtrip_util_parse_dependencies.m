function test_suite=test_moxunit_fieldtrip_util_parse_dependencies
% #   For MOxUnit_fieldtrip's copyright information and license terms,   #
% #   see the COPYING file distributed with MOxUnit_fieldtrip.           #

    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function varargout=randint(mx)
    n=nargout;
    varargout=arrayfun(@(unused)ceil(rand()*mx),1:n,'UniformOutput',false);

function assert_output(input, expected_output)
    output=moxunit_fieldtrip_util_parse_dependency(input);
    if isempty(output)
        assertTrue(isempty(expected_output));
    else
        assertEqual(output,expected_output,sprintf('input: ''%s''',input));
    end

function test_fieldtrip_util_parse_dependencies_basics()
    empty={};
    input_outputs={'foo',{'foo'};...
                  'foo  ',{'foo'};...               % whitespace
                  'foo bar',{'foo','bar'};...
                  'foo  bar  baz ',{'foo','bar','baz'};...
                  '',empty;...
                  };
    for k=1:size(input_outputs)
        input_output=input_outputs(k,:);
        assert_output(input_output{1},input_output{2});
    end

function moxunit_fieldtrip_util_parse_dependency_exceptions
    aet=@(varargin)assertExceptionThrown(@()...
            moxunit_fieldtrip_util_parse_dependency(varargin{:}),'');
    % inputs that are not string throw an exception
    aet(struct)
    aet([1 2])
    aet(1);
