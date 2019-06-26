function test_suite=test_moxunit_fieldtrip_util_struct2json
% #   For MOxUnit_fieldtrip's copyright information and license terms,   #
% #   see the COPYING file distributed with MOxUnit_fieldtrip.           #

    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function s=randstr()
    s=char(ceil(rand(1,10)*26+64));


function test_moxunit_fieldtrip_util_struct2json_basics
    s=randstr();
    t=randstr();
    d=rand();
    arg_output={{},'{}';...
                {s,d},sprintf('{"%s":%d}',s,d);...
                {s,t},sprintf('{"%s":"%s"}',s,t);...
                {s,t,t,d},sprintf('{"%s":"%s","%s":%d}',s,t,t,d);...
                {s,sprintf('"%s"X"',t)},...
                                sprintf('{"%s":"\\"%s\\"X\\""}',s,t);...
                {s,true},sprintf('{"%s":true}',s);...
                {s,false},sprintf('{"%s":false}',s);...

                };

    n=size(arg_output,1);
    for k=1:n
        arg=arg_output{k,1};

        input=struct(arg{:});
        output=arg_output{k,2};

        assertEqual(moxunit_fieldtrip_util_struct2json(input),output);
    end





function test_moxunit_fieldtrip_util_struct2json_exceptions
    aet=@(varargin)assertExceptionThrown(@()...
                moxunit_fieldtrip_util_struct2json(varargin{:}),'');

    bad_inputs={'foo',...
                3,...
                struct('k',[1 2])};
    for k=1:numel(bad_inputs)
        aet(bad_inputs{k});
    end
