function test_suite=test_moxunit_fieldtrip_test_suite_filtering
% #   For MOxUnit_fieldtrip's copyright information and license terms,   #
% #   see the COPYING file distributed with MOxUnit_fieldtrip.           #

    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;


function test_fieldtrip_suite_filter_basics()
% test that object can be instantiated and has expected properties

    % each element in the struct below containes these elements:
    % - the label for suite.filter ('maxwalltime','maxmem','loadfile')
    % - the value for suite.filter
    % - lines for which test cases should be kept    after the filter
    % - "                                  " removed "              "
    s=struct();
    s.walltime={'maxwalltime',...
                {20,'20'},...
                {'','% WALLTIME 10','% WALLTIME 00:00:10'},...
                {'% WALLTIME 30','% WALLTIME 00:05:00'}};

    s.mem={'maxmem',...
                {10000,'10000'},...
                {'','% MEM 10','% MEM 1kb'},...
                {'% MEM 20000','% MEM 1tb'}};
    s.dccnpath={'loadfile',...
                {false,'no'},...
                {'','% dccnpath'},...
                {'dccnpath','fullfile(dccnpath,''bar'')'}};

    keys=fieldnames(s);
    for k=1:numel(keys)
        key=keys{k};
        prop=s.(key);

        label=prop{1};
        values=prop{2};

        for skip_test=[false,true]
            prop_idx=3;
            if skip_test
                prop_idx=prop_idx+1;
            end
            line_cell=prop{prop_idx};

            for j=1:numel(line_cell)
                line=line_cell{j};
                for m=1:numel(values)
                    value=values{m};

                    test_dir=tempname();
                    fn=randtestfilename(test_dir);

                    helper_test_single_test_case_filter(fn,...
                                                        label,...
                                                        value,...
                                                        line,...
                                                        skip_test);
                end
            end
        end
    end

function test_fieldtrip_suite_filter_exclude_if_prefix_equals_failed()
    test_dir=tempname();

    prefix='failed_';

    arg='exclude_if_prefix_equals_failed';
    no_yes_values={{'no',false},...
                    {'yes',true}};
    for k=1:numel(no_yes_values)
        exclude=k>1;
        for j=1:numel(no_yes_values{k})
            value=no_yes_values{k}{j};

            fn=randtestfilename(test_dir,prefix);
            helper_test_single_test_case_filter(fn,arg,value,'',exclude);
        end
    end



function helper_test_single_test_case_filter(fn,arg,value,...
                                                line,...
                                                skip_test)

    % make test case
    test_dir=fileparts(fn);
    suite=make_basic_test_suite(test_dir);
    dir_cleaner=onCleanup(@()clean_dir(test_dir));

    tc=build_test_case(fn, {line});

    % suite is empty initially
    assertEqual(countTestCases(suite),0);

    % after adding the test it must be in the suite
    suite=addFromFile(suite,fn);
    assertEqual(countTestCases(suite),1);
    assertEqual(getLocation(getTestNode(suite,1)),...
                    getLocation(tc));

    % skip test
    suite=filter(suite,arg,value);

    % tests must still be there
    assertEqual(countTestCases(suite),1);

    report=MOxUnitTestReport(0,1,'');
    report=run(suite,report);

    assertEqual(countTestOutcomes(report),1);

    % get test outcome.
    outcome=getTestOutcome(report,1);

    % test should not have failed
    assertEqual(isNonFailure(outcome),true);

    % succesful if test was not skipped
    is_success=isSuccess(outcome);
    assertEqual(is_success,~skip_test);

    if skip_test
        assertTrue(isa(outcome,'MOxUnitSkippedTestOutcome'));
    end



function s=randstr
    s=char(ceil(rand(1,10)*26+64));


function suite=make_basic_test_suite(test_dir)
    suite=MOxUnitFieldTripTestSuite();
    mkdir(test_dir);


function fn=randtestfilename(parent_dir,prefix)
    if nargin<2
        prefix='test_';
    end
    test_name=randstr();
    func_name=sprintf('%s%s',prefix,test_name);
    fn=fullfile(parent_dir,sprintf('%s.m',func_name));


function clean_dir(dir_name)
    ds=dir(dir_name);
    n=numel(ds);
    for k=1:n
        d=ds(k);
        if ~d.isdir
            % delete file
            pth=fullfile(dir_name,d.name);
            delete(pth);
        end
    end

    rmdir(dir_name);


function tc=build_test_case(fn, body_lines)
    [unused,func_name]=fileparts(fn);
    lines=[sprintf('function %s',func_name),...
            body_lines];
    write_content(fn,lines);
    tc=MOxUnitFieldTripTestCase(fn);


function write_content(fn,lines)
    fid=fopen(fn,'w');
    cleaner=onCleanup(@()fclose(fid));
    fprintf(fid,'%s\n',lines{:});
