%TEST_ALL Script for running most of the demos in GPstuff and comparing new values
%         with previously saved ones. 
%
%  Description
%    Run each test function in xunit folder. Each of the test functions
%    runs corresponding demo and compares the values to previously saved
%    ones for errors. This is useful e.g. in case user modifies functions
%    provided by GPstuff.
%
%  See also
%    Readme.txt

% Copyright (c) 2011-2012 Ville Tolvanen

% Remove previous test values and create empty folder for new test values.
path = which('test_all.m');
path = strrep(path, 'test_all.m', 'testValues');
if exist(path,'dir') == 7
    rmdir(path, 's');
end
mkdir(path);

% Create TestSuite object suite which includes all tests from root
% directory and logger object which keeps track of the number of tests, failures and
% errors. Run tests and print logger information. Can be run from any
% directory.

path2 = cd;
cd(strrep(path,'testValues', ''));
suite = TestSuite.fromPwd();

logger = TestRunLogger;
suite.run(logger);
cd(path2);

logger
