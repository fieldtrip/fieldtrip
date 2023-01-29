function errors = testStringServices(host)
%% Shows how to call hed-services to process a list of hedstrings.
% 
%  Example 1: Validate valid list of strings using HED version.
%
%  Example 2: Validate invalid list of strings using HED URL.
%
%  Example 3: Validate invalid list of strings uploading HED schema.
%
%  Example 4: Convert valid strings to long using HED version.
%

%% Get the options and data
[servicesUrl, options] = getHostOptions(host);
data = getTestData();
errors = {};


%% Example 1: Validate valid list of strings using HED URL.
request1 = struct('service', 'strings_validate', ...
                  'schema_version', '8.0.0', ...
                  'string_list', '', ...
                  'check_warnings_validate', 'on');
request1.string_list = data.goodStrings;
response1 = webwrite(servicesUrl, request1, options);
response1 = jsondecode(response1);
outputReport(response1, 'Example 1 Validating a valid list of strings');
if ~isempty(response1.error_type) || ...
   ~strcmpi(response1.results.msg_category, 'success')
   errors{end + 1} = 'Example 1 failed to validate valid HED strings.';
end

%% Example 2: Validate a list of invalid strings. HED schema is URL.
request2 = struct('service', 'strings_validate', ...
                  'schema_url', '', ...
                  'string_list', '', ...
                  'check_for_warnings', 'on');
request2.string_list = data.badStrings;
request2.schema_url = data.schemaUrl;
response2 = webwrite(servicesUrl, request2, options);
response2 = jsondecode(response2);
outputReport(response2, ...
             'Example 2 validating a list of strings with invalid values');
if ~isempty(response2.error_type) || ...
   ~strcmpi(response2.results.msg_category, 'warning')
   errors{end + 1} = 'Example 2 failed to detect invalid HED strings.';
end

%% Example 3: Validate list of invalid strings uploading HED schema.
request3 = struct('service', 'strings_validate', ...
                  'schema_string', data.schemaText, ...
                  'string_list', '', ...
                  'check_for_warnings', 'on');
request3.string_list = data.badStrings;               
response3 = webwrite(servicesUrl, request3, options);
response3 = jsondecode(response3);
outputReport(response3, ...
            'Example 3 validate invalid strings using an uploaded HED schema');
if ~isempty(response3.error_type) || ...
   ~strcmpi(response3.results.msg_category, 'warning')
   errors{end + 1} = 'Example 3 failed to detect invalid HED strings.';
end

%% Example 4: Convert valid strings to long using HED version.
request4 = struct('service', 'strings_to_long', ...
                  'schema_version', '8.0.0', ...
                  'string_list', '');
request4.string_list = data.goodStrings;               
response4 = webwrite(servicesUrl, request4, options);
response4 = jsondecode(response4);
outputReport(response4, ...
             'Example 4 Convert a list of valid strings to long');
if ~isempty(response4.error_type) || ...
   ~strcmpi(response4.results.msg_category, 'success')
   errors{end + 1} = 'Example 4 failed to convert HED strings for long.';
end

