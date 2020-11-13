function [SortedStruct, SortedIndeces] = sortNEV(STRUCT, field, TestValue, Report)

% This function extracts information out of a structure based on the field
% values. The outputs can be the sorted structure and the indeces of those
% extracted fields. 
%
% Use [SortedStruct, SortedIndeces] = sortNEV(STRUCT, field, TestValue, Report)
%
%   STRUCT:     The structure to be used for extraction.
%
%   field:      The field in STRUCT that is set as criteria for extraction.
%               This parameter is optional. If not passed, the function
%               will prompt the user to provide the name of the field.
%
%   TestValue:  The value that "field" needs to be equal to in order to
%               get extracted.
%               This parameter is optional. If not passed, the function
%               will prompt the user to provide the value.
%
%   Report:     If this flag is set to 1 the function will show a short
%               summary of the data that was processed.
%               This parametere is optional.
%               DEFAULT: will not show report.
%
%   OUTPUT
%
%   SortedStruct:   The sorted structure output.
%   
%   SortedIndeces:  The indeces to the sorted elements.
%
%   Example: 
%   
%   [NewStruct, Indeces] = sortNEV(MainStruct, 'TestField', 3, 1);
%
%   In the example above, the MainStruct is the structure containing all
%   elements. The function will search through MainStruct and extract all
%   elements that have their 'TestField' field set to '3'. It will also
%   show a report of how many elemenets were processed and how many were
%   extracted.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   ktorab@blackrockmicro.com
%   Blackrock Microsystems
%   
%   Version 1.4.1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Validates all passed on parameters and will prompt user for any missing
%  parameter.

if ~exist('field', 'var')
    field = input('What is the name of the field of interest? ');
end

while ~isfield(STRUCT, field)
    display('The field name does not exist. Try again.');
    field = input('What is the name of the field of interest? ');
end

if ~exist('TestValue', 'var')
    TestValue = input('What value does the field need to be equal to (string format)? ');
end

if ~exist('Report', 'var')
    Report = 0;
end

%% Searches through the structure to extract the elements of interest.
SortedIndeces = NaN(1,length(STRUCT));
if isa(TestValue, 'double')
%     display(['Please verify to make sure the field ' field ' does not contain any [] data or the sort will not be successful.'])
    for strIDX = 1:length(STRUCT)
        if find(STRUCT(strIDX).(field) == TestValue)
            SortedIndeces(strIDX) = 1;
        else
            SortedIndeces(strIDX) = 0;
        end
    end
    SortedStruct = STRUCT(find(SortedIndeces==1));
else
    for strIDX = 1:length(STRUCT)
        if strcmpi(STRUCT(strIDX).(field), TestValue)
            SortedIndeces(strIDX) = 1;
        else
            SortedIndeces(strIDX) = 0;
        end
    end
    SortedStruct = STRUCT(find(SortedIndeces==1));
end


%% Will define as an empty structure of no instances were found.
if ~exist('SortedStruct', 'var')
    display('No instances found.');
    SortedStruct = struct([]);
else

%% If the flag is set, it will show a report of how many files were
%  processed.
if Report == 1
    if isempty(SortedStruct)
        display('No instances found.');
    else
        display([num2str(length(SortedStruct)) ' many instances were found and extracted.']);
    end
end

%% Will force the function to send an output argument to the caller
%  environment.
if (nargout == 0)
    assignin('caller', 'SortedStruct', SortedStruct);
end

end