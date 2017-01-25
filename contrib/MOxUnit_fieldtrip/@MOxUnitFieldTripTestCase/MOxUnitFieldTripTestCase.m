function obj=MOxUnitFieldTripTestCase(fn)
% Initialize MOxUnitFieldTripTestCase instance
%
% obj=MOxUnitFieldTripTestCase(filename)
%
% Inputs:
%   filename                Filename of .m file with FieldTrip test
%
% Output:
%   obj                     MOxUnitFieldTripTestCase representing the
%                           case in the filename.
%
% Notes:
%  - unlike its MOxUnitFunctionHandleTestCase superclass, this class
%    represents a test case through a single function in an mfile.
%    Running the tests involves 'cd'-ing to the directory where the mfile
%    is stored, and then calling the function.
%  - the 'get' function can be used to get 'mem', 'walltime', or 'dccnpath'
%    properties
%
% #   For MOxUnit_fieldtrip's copyright information and license terms,   #
% #   see the COPYING file distributed with MOxUnit_fieldtrip.           #

    full_path=make_full_path(fn);
    if ~exist(full_path,'file')
        error('file ''%s'' does not exist',fn);
    end

    s=struct();

    % get 'walltime', 'mem' and 'dccnpath' properties, if present
    s.properties=moxunit_fieldtrip_util_read_test_mfile_properties(fn);

    [parent_dir,name]=fileparts(full_path);
    function_handle=@()run_function_in_dir(parent_dir,name);
    location=full_path;

    obj=class(s,'MOxUnitFieldTripTestCase',...
                    MOxUnitFunctionHandleTestCase(...
                                        name,location,function_handle));


function run_function_in_dir(in_dir,name)
    % This is a wrapper function used to run each FieldTrip test

    % If the test changes the path, make sure to undo such changes
    % afterwards
    orig_path=path();
    path_restorer=onCleanup(@()path(orig_path));

    % Some tests open new figures and do not close them afterwards, which
    % can lead to a proliferation of figure windows when running a suite
    % with many tests. Here we want to close any figures that were opened
    % while running the test.
    % The current approach uses 'findobj' to find the figures open both
    % before and after running the test, and then closes the ones that were
    % present in the latter but not the former list using 'setdiff'.
    % This approach works under Octave and Matlab, even though findobj
    % returns a numeric vector in Octave and a Figure object in Matlab.
    orig_fig_handles=findobj('Type','figure');
    fig_cleaner=onCleanup(@()close(setdiff(findobj('Type','figure'),...
                                                    orig_fig_handles)));

    % keep track of the original working directory, and ensure that we
    % return to that directory after that test
    orig_dir=pwd();
    pwd_restorer=onCleanup(@()cd(orig_dir));

    % run the test
    cd(in_dir);
    handle=str2func(name);
    handle();

    % flush event queue
    drawnow();


function full_path=make_full_path(fn)
    if is_absolute_path(fn)
        prefix='';
    else
        prefix=pwd();
    end

    full_path=fullfile(prefix,fn);

function tf=is_absolute_path(fn)
    if ispc()
        % second character must be ':'
        abs_path_pattern='^[a-zA-Z]:';
    else
        % first character must be '/'
        abs_path_pattern='^/';
    end

    tf=~isempty(regexp(fn,abs_path_pattern,'once'));
