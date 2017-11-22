function test_suite=test_moxunit_fieldtrip_test_case
% #   For MOxUnit_fieldtrip's copyright information and license terms,   #
% #   see the COPYING file distributed with MOxUnit_fieldtrip.           #

    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;


function test_fieldtrip_test_case_basics
% test that object can be instantiated and has expected properties
    fn=randtestfilename();
    [unused,func_name]=fileparts(fn);

    % before writing lines should give an error
    constructor=@()MOxUnitFieldTripTestCase(fn);
    assertExceptionThrown(constructor,'');

    % write contents
    lines={sprintf('function %s',func_name),...
                        '% do nothing'};
    write_content(fn,lines);
    cleaner=onCleanup(@()delete(fn));

    tc=constructor();

    assert(isa(tc,'MOxUnitTestCase'));
    assertEqual(getName(tc),func_name);
    assertEqual(getLocation(tc),fn);


function test_fieldtrip_test_case_success()
% simple passing test
    helper_test_case_with_content('abs(2);',true);


function test_fieldtrip_test_case_failure()
% simple failing test
    helper_test_case_with_content('error(''here'');',false);


function test_fieldtrip_test_case_error()
% simple erroring test
    helper_test_case_with_content('wrong))',false);


function test_moxunit_fieldtrip_test_suite_basics()
% simple instantiation
    suite=MOxUnitFieldTripTestSuite();
    assertTrue(isa(suite,'MOxUnitTestSuite'));
    assertEqual(countTestCases(suite),0);


function test_moxunit_fieldtrip_test_suite_success()
% simple passing test
    helper_test_suite_with_content('abs(2);',true);


function test_fieldtrip_test_suite_failure()
% simple failing test
    helper_test_suite_with_content('error(''here'');',false);


function test_fieldtrip_test_suite_error()
% simple erroring test
    helper_test_suite_with_content('wrong))',false);

%%%%%%%%%%%%%%%%%%
% helper functions


function s=randstr
    s=char(ceil(rand(1,10)*26+64));

function fn=randtestfilename(parent_dir)
    if nargin<1
        parent_dir=tempdir();
    end
    test_name=randstr();
    func_name=sprintf('test_%s',test_name);
    fn=fullfile(parent_dir,sprintf('%s.m',func_name));



function helper_test_case_with_content(single_line,should_pass)
    fn=randtestfilename();

    cleaner=onCleanup(@()delete(fn));
    tc=helper_build_test_case(fn, single_line);

    report=MOxUnitTestReport(0,1,'');
    report=run(tc,report);

    result=wasSuccessful(report);
    assertEqual(result,should_pass);


function helper_test_suite_with_content(single_line,should_pass)
    test_dir=tempname();
    mkdir(test_dir);
    dir_cleaner=onCleanup(@()helper_clean_dir(test_dir));

    fn=randtestfilename(test_dir);
    tc=helper_build_test_case(fn, single_line);

    % add tests either from file or from test object
    test_adders={@(s) addFromFile(s,fn),...
                 @(s) addTest(s,tc)};

    for k=1:numel(test_adders)
        adder=test_adders{k};

        suite=MOxUnitFieldTripTestSuite();
        suite=adder(suite);
        assertEqual(countTestCases(suite),1);

        report=MOxUnitTestReport(0,1,'');
        report=run(suite,report);

        result=wasSuccessful(report);
        assertEqual(result,should_pass);
    end


function helper_clean_dir(dir_name)
    ds=dir(dir_name);
    n=numel(ds);
    for k=1:n
        d=ds(k);
        if ~d.isdir
            pth=fullfile(dir_name,d.name);
            delete(pth);
        end
    end

    rmdir(dir_name);


function tc=helper_build_test_case(fn, single_line)
    [unused,func_name]=fileparts(fn);
    lines={sprintf('function %s',func_name),...
            single_line};
    write_content(fn,lines);
    tc=MOxUnitFieldTripTestCase(fn);


function write_content(fn,lines)
    fid=fopen(fn,'w');
    cleaner=onCleanup(@()fclose(fid));
    fprintf(fid,'%s\n',lines{:});
