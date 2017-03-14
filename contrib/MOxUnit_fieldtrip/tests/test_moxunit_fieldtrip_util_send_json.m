function test_suite=test_moxunit_fieldtrip_util_send_json
% #   For MOxUnit_fieldtrip's copyright information and license terms,   #
% #   see the COPYING file distributed with MOxUnit_fieldtrip.           #

    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function s=randstr()
    s=char(ceil(rand(1,10)*26+64));


function test_moxunit_fieldtrip_util_send_json_methods
    methods=moxunit_fieldtrip_util_send_json();
    assert(iscellstr(methods));

    assertEqual(isempty(which('webwrite')),...
                        all(isempty(strmatch('webwrite',methods))));
    if ispc()
        % assume that curl is not present
        assertTrue(all(isempty(strmatch('curl',methods))));
    else
        assertEqual(isempty(which('curl')),...
                            unix('which curl>/dev/null')~=0);
    end


function test_moxunit_fieldtrip_util_struct2json_exceptions
    aet=@(varargin)assertExceptionThrown(@()...
                moxunit_fieldtrip_util_struct2json(varargin{:}),'');


    bad_inputs={{'foo'},...
                 {3,struct},...
                {'foo','bar'}};
    for k=1:numel(bad_inputs)
        args=bad_inputs{k};
        aet(args);
    end
