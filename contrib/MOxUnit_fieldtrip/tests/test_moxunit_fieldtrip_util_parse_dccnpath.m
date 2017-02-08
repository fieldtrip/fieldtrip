function test_suite=test_moxunit_fieldtrip_util_parse_dccnpath
% #   For MOxUnit_fieldtrip's copyright information and license terms,   #
% #   see the COPYING file distributed with MOxUnit_fieldtrip.           #

    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;


function assert_output(input, expected_output)
    output=moxunit_fieldtrip_util_parse_dccnpath(input);
    assertEqual(output,expected_output);

function test_fieldtrip_util_parse_dccnpath_basics
    % just the statement itself
    assert_output('dccnpath',true);

    % with whitespace
    assert_output('   dccnpath   (',true);

    % not if in a comment
    assert_output('   % dccnpath   (',[]);

    % not in a string
    assert_output('   foo(''dccnpath'')',[]);

    % not part of another string or in parts
    assert_output('   dccnpathsuffix',[]);
    assert_output('   prefixdccnpath',[]);
    assert_output('   _dccnpath ',[]);
    assert_output('   dccnpath_ ',[]);
    assert_output('   dccnpath1 ',[]);
    assert_output('   dccn path   ',[]);




function test_fieldtrip_util_parse_walltime_exceptions
    aet=@(varargin)assertExceptionThrown(@()...
            moxunit_fieldtrip_util_parse_dccnpath(varargin{:}),'');
    % inputs that are not string throw an exception
    aet(struct)
    aet([1 2])
    aet(1);
