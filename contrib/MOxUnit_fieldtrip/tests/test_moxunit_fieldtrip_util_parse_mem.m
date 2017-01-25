function test_suite=test_moxunit_fieldtrip_util_parse_mem
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
    output=moxunit_fieldtrip_util_parse_mem(input);
    assertEqual(output,expected_output);

function test_fieldtrip_util_parse_walltime_str
    suffixes={'kb','mb','gb','tb',''};
    for k=1:numel(suffixes)
        suffix=suffixes{k};
        v=randint(10000);

        input=sprintf('%d%s',v,suffix);

        if isempty(suffix)
            mply=1;
        else
            mply=2^(k*10);
        end

        assert_output(input,mply*v);

        % not in a comment, not valid
        input=sprintf('% %d%s',v,suffix);
        assert_output(input,[]);
    end

    % test infinity
    assert_output('inf',inf);

    % test more whitespace
    assert_output(' 3 % comment',3);


function test_fieldtrip_util_parse_walltime_exceptions
    aet=@(varargin)assertExceptionThrown(@()...
            moxunit_fieldtrip_util_parse_mem(varargin{:}),'');
    % inputs that are not string throw an exception
    aet(struct)
    aet([1 2])
    aet(1);
