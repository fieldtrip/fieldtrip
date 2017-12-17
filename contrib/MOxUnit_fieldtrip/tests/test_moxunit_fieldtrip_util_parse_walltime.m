function test_suite=test_moxunit_fieldtrip_util_parse_walltime
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
    output=moxunit_fieldtrip_util_parse_walltime(input);
    assertEqual(output,expected_output);

function test_fieldtrip_util_parse_walltime
    % seconds
    many_ss=randint(1000);
    assert_output(sprintf('%d',many_ss),...
                                many_ss);

    % hours, minutes, seconds
    [hh,mm,ss]=randint(60);
    assert_output(sprintf('%02d:%02d:%02d',hh,mm,ss),...
                                (hh*60+mm)*60+ss);

    % not a number
    [hh,mm,ss]=randint(60);
    assert_output('foo',[]);

    % should deal with comments after
    assert_output('  3 % comment',3);



function test_fieldtrip_util_parse_walltime_exceptions
    aet=@(varargin)assertExceptionThrown(@()...
            moxunit_fieldtrip_util_parse_walltime(varargin{:}),'');
    % inputs that are not string throw an exception
    aet(struct)
    aet([1 2])
    aet(1);


