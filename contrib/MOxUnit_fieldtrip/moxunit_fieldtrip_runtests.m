function passed=moxunit_fieldtrip_runtests(varargin)
% run Fieldtrip tests with MOxUnit framework
%
% passed=moxunit_fieldtrip_runtests(...)
%
% Inputs:
%   'maxwalltime', w        Maximum duration as may be specified in each
%                           test; either numeric in seconds, or a string
%                           in format HH:MM:SS
%   'maxmem',m              Maximum memory usage as may be specified in
%                           each test; either numeric in bytes, or a string
%                           with a number and suffix 'kb', 'mb', 'gb', or
%                           'tb'
%   'upload',u              If 'true' or 'yes', then results are uploaded
%                           to the FieldTrip dashboard (default: true if
%                           the platform provides the 'weboptions' and
%                           'weboptions' functions)
%   'dependency',d          Do not run test if it matches a dependency in d
%   'xmloutput',f           Write JUnit-like XML file with test results to
%                           file f. This can be used with the shippable.com
%                           continuous integration service.
%   'exclude_if_prefix_equals_failed',e
%                           Exclude the test if the filename starts with
%                           'failed'.

%
% Output:
%   passed                  scalar boolean that is false if one or more
%                           tests failed or errored, and true otherwise.
%                           Tests that were skipped
%
% Examples:
%   % run all tests
%   moxunit_fieldtrip_runtests
%   %
%   % run tests with an expected execuation time of at most 10 minutes,
%   % maxmimum memory usage of 1 gigabyte, and with upload disabled.
%   moxunit_fieldtrip_runtests maxwalltime 00:10:00 maxmem 1gb upload no
%   %
%   % run all tests and store JUnit test results in file foo.xml. This file
%   % contains for each tests whether it passed, failed, errored, or was
%   % skipped; if it failed or errored, the full stack trace; if it was
%   % skipped, the reason why it was skipped.
%   moxunit_fieldtrip_runtests xmloutput foo.xml
%
% Notes:
% - This function aims to provide similar syntax as 'ft_test run', but
% provides an output
%
% #   For MOxUnit_fieldtrip's copyright information and license terms,   #
% #   see the COPYING file distributed with MOxUnit_fieldtrip.           #
%
% See also: ft_test

    passed=false; % fail unless all tests were successful

    % parse arguments
    opt=parse_args(varargin{:});

    % make a suite and add the tests
    suite=MOxUnitFieldTripTestSuite();
    suite=add_test_files(suite,opt.filelist);

    % show number of tests
    disp(suite);

    % apply filter args
    suite=filter(suite,opt.filter_args{:});

    % set output options
    verbosity=2;     % show detailed information
    output_stream=1; % standard output

    n_tests=countTestCases(suite);
    if n_tests==0
        fprintf(output_stream,['No test cases found - this is '...
                            'considered as a test suite failure\n']);
        return;
    end



    report=MOxUnitFieldTripTestReport(verbosity,output_stream);
    report=run(suite,report);

    if opt.upload
        sendToFieldtripDashboard(report);
    else
        fprintf(output_stream,'Not sending results to dashboard\n');
    end

    if ~isempty(opt.xmloutput)
        writeXML(report,opt.xmloutput);
        fprintf(output_stream,'XML test result written to ''%s''\n',...
                                                opt.xmloutput);
    end

    disp(report);
    passed=wasSuccessful(report);


function suite=add_test_files(suite,filelist)
% add tests specifed in filelist. If no tests are specified, then all
% tests in the 'test' directory are added
    if isempty(filelist)
        fieldtrip_test_dir=get_fieldtrip_test_dir();
        fieldtrip_test_file_pattern='^(test|failed).*\.m$';
        suite=addFromDirectory(suite,fieldtrip_test_dir,...
                                    fieldtrip_test_file_pattern);
    else
        % add tests from the FieldTrip test directory, so that file
        % names of tests do not need the full path to be added
        orig_pwd=pwd();
        pwd_resetter=onCleanup(@()cd(orig_pwd));
        cd(get_fieldtrip_test_dir());

        for k=1:numel(filelist)
            suite=addFromFile(suite,filelist{k});
        end
    end


function opt=parse_args(varargin)
    filter_args={};
    upload=[];
    filelist={};
    loadfile=[];
    xmloutput=[];

    pos=0;
    while true
        pos=pos+1;

        if pos>nargin
            break
        end

        key=varargin{pos};

        switch key
            case 'upload'
                [pos,value]=get_next_element(varargin,pos);
                upload=istrue(value);

            case 'loadfile'
                [pos,value]=get_next_element(varargin,pos);
                loadfile=istrue(value);

            case {'dependency','maxmem','maxwalltime',...
                        'exclude_if_prefix_equals_failed'}
                [pos,value]=get_next_element(varargin,pos);
                filter_args{end+1}=key;
                filter_args{end+1}=value;

            case 'xmloutput'
                [pos,value]=get_next_element(varargin,pos);
                xmloutput=value;

            otherwise
                filelist{end+1}=key;

        end
    end

    if isempty(loadfile)
        % 'loadfile' was not set as argument
        % set to true if at DCCN-like computer, false otherwise
        loadfile=exist(get_fieldtrip_data_dir(),'dir');
    end

    if isempty(upload)
        % TODO: decide whether upload should be disabled when not running
        % from DCCN-like computer

        % See if there are any upload methods available
        upload=~isempty(moxunit_fieldtrip_util_send_json());
    end


    filter_args=[filter_args,{'loadfile',loadfile}];

    opt=struct();
    opt.filter_args=filter_args;
    opt.upload=upload;
    opt.filelist=filelist;
    opt.xmloutput=xmloutput;


function [pos,value]=get_next_element(argcell,pos)

    if pos+1>numel(argcell)
        error('Missing value after key ''%s''', argcell{pos});
    end

    pos=pos+1;
    value=argcell{pos};

function result=ft_platform_supports_wrapper(arg)
    result=fieldtrip_function_wrapper(fullfile('utilities','private'),...
                                        'ft_platform_supports',arg);

function result=fieldtrip_function_wrapper(subdir,func_name,varargin)
    orig_pwd=pwd();
    cleaner=onCleanup(@()cd(orig_pwd));

    eval_dir=fullfile(get_fieldtrip_root_dir(),subdir);
    cd(eval_dir);

    func=str2func(func_name);
    result=func(varargin{:});



function fieldtrip_root_dir=get_fieldtrip_root_dir()
    fieldtrip_root_dir=fileparts(which('ft_defaults'));
    if isempty(fieldtrip_root_dir)
        error('Fieldtrip path is not set');
    end

function fieldtrip_test_dir=get_fieldtrip_test_dir()
    fieldtrip_root_dir=get_fieldtrip_root_dir();
    fieldtrip_test_dir=fullfile(fieldtrip_root_dir,'test');

function fieldtrip_data_dir=get_fieldtrip_data_dir()
    fieldtrip_data_dir='/home/common/matlab/fieldtrip/data/';
    fieldtrip_data_dir=fieldtrip_function_wrapper('test',...
                                'dccnpath',fieldtrip_data_dir);


