function passed = ft_test_moxunit_run(unused,varargin)

% FT_TEST_MOXUNIT_RUN

% Copyright (C) 2017, Nikolaas N. Oosterhof
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

    % ensure path is set for MOxUnit fieldtrip test functions, but
    % set to original state after running this function
    orig_path=path();
    path_resetter=onCleanup(@()path(orig_path));
    ensure_moxunit_path_is_set();
    check_dependencies();

    % parse arguments
    [filter_args,upload,filelist]=parse_args(varargin{:});

    % make a suite and add the tests
    suite=MOxUnitFieldTripTestSuite();
    suite=add_test_files(suite,filelist);

    % show number of tests
    disp(suite);

    % apply filter args
    suite=filter(suite,filter_args{:});

    verbosity=2;
    output_stream=1;
    report=MOxUnitFieldTripTestReport(verbosity,output_stream);
    report=run(suite,report);

    if upload
        send_to_fieldtrip_dashboard(report);
    else
        fprintf(output_stream,'Not sending results to dashboard\n');
    end

    disp(report);
    passed=wasSuccessful(report);


function suite=add_test_files(suite,filelist)
% add tests specifed in filelist. If no tests are specified, then all
% tests in the 'test' directory are added
    if isempty(filelist)
        fieldtrip_test_dir=get_fieldtrip_test_dir();
        fieldtrip_test_file_pattern='^[Tt][Ee][sS][tT].*\.m$';
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


function [filter_args,upload,filelist]=parse_args(varargin)
    filter_args={};
    upload=[];
    filelist={};
    loadfile=[];

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

            case {'dependency','maxmem','maxwalltime'}
                [pos,value]=get_next_element(varargin,pos);
                filter_args{end+1}=key;
                filter_args{end+1}=value;


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
        % TODO: decide whehter upload should be disabled when not running
        % from DCCN-like computer
        upload=ft_platform_supports_wrapper('weboptions') && ...
                ft_platform_supports_wrapper('webwrite');
    end


    filter_args=[filter_args,{'loadfile',loadfile}];


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


function check_dependencies()
% throw an error if MOxUnit is not available

    if isempty(which('moxunit_runtests'))
        error(['MOxUnit is required; see '...
                    'https://github.com/moxunit/moxunit']);
    end

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


function ensure_moxunit_path_is_set()
    if isempty(which('MOxUnitFieldTripTestSuite'))
        moxunit_fieldtrip_dir=fullfile(get_fieldtrip_root_dir(),...
                                       'contrib','MOxUnit_fieldtrip');
        addpath(moxunit_fieldtrip_dir);
    end
