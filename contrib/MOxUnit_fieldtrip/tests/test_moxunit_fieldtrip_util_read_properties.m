function test_suite=test_moxunit_fieldtrip_util_read_properties
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

function s=randstr
    s=char(ceil(rand(1,10)*26+64));


function test_fieldtrip_test_case_with_walltime
    many_ss=randint(1e4);
    helper_test_property(sprintf('%% WALLTIME %d',many_ss),...
                        'walltime',many_ss);


    [hh,mm,ss]=randint(60);
    helper_test_property(sprintf('%% WALLTIME %02d:%02d:%02d',hh,mm,ss),...
                        'walltime',(hh*60+mm)*60+ss);

    % if not in a comment then it is not set
    helper_test_property(sprintf('WALLTIME %02d:%02d:%02d',hh,mm,ss),...
                        'walltime',[]);


function test_fieldtrip_test_case_with_memory
    [tb,gb,mb,kb,b]=randint(100);
    helper_test_property(sprintf('%% MEM %dtb',tb),'mem',tb*2^40);
    helper_test_property(sprintf('%% MEM %dgb',gb),'mem',gb*2^30);
    helper_test_property(sprintf('%% MEM %dmb',mb),'mem',mb*2^20);
    helper_test_property(sprintf('%% MEM %dkb',kb),'mem',kb*2^10);
    helper_test_property(sprintf('%% MEM %d',b),  'mem', b);

    % if not in a comment then it is not set
    helper_test_property(sprintf('MEM %d',b),  'mem', []);

function test_fieldtrip_test_case_with_dccnpath
    helper_test_property('','dccnpath',[]);
    helper_test_property('dccnpath % unused','dccnpath',true);

function test_fieldtrip_test_case_with_dependency
    s=randstr();
    t=randstr();

    helper_test_property(sprintf('%% TEST %s',s),'dependency',{s});
    helper_test_property(sprintf('%% TEST %s %s',s,t),'dependency',{s,t});
    helper_test_property(sprintf('%% TEST %s %s ',s,t),'dependency',{s,t});

    % if not in a comment then it is not set
    helper_test_property(sprintf('TEST %s',s),'dependency',[]);


function test_fieldtrip_test_case_no_properties
    keys={'mem','walltime','dccnpath','dependency','foo'};

    fn=helper_write_test_case('unused');
    cleaner=onCleanup(@()delete(fn));

    % test both the moxunit_fieldtrip_util_read_mfile_properties(fn)
    % function and the constructor for test cases
    s=moxunit_fieldtrip_util_read_test_mfile_properties(fn);
    tc=MOxUnitFieldTripTestCase(fn);

    % when property is not set, it should not be present in the test case
    for k=1:numel(keys)
        key=keys{k};
        assert(~isfield(s,key));
        assert(isempty(get(tc,key)))
    end


function fn=helper_write_test_case(content)
    fn=sprintf('%s.m',tempname());
    fid=fopen(fn,'w');
    cleaner=onCleanup(@()fclose(fid));

    lines={'function test_me',...
            sprintf('%s',content),...
            'abs(2)'};

    fprintf(fid,'%s\n',lines{:});


function helper_test_property(content, key, expected_value)
    fn=helper_write_test_case(content);
    cleaner=onCleanup(@()delete(fn));

    s=moxunit_fieldtrip_util_read_test_mfile_properties(fn);
    if isempty(expected_value)
        assert(~isfield(s,key));
    else
        value=s.(key);
        assertEqual(value, expected_value);
    end

    tc=MOxUnitFieldTripTestCase(fn);
    value=get(tc,key);
    assertEqual(value, expected_value);




