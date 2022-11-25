function errors = testEventServices(host)
%% Shows how to call hed-services to searcg a BIDS events file.
%
%  Example 1: Search an events file for HED using a valid query.
%
%  Example 2: Search an events file for HED and return additional columns.

%% Get the options and data
[servicesUrl, options] = getHostOptions(host);
data = getTestData();
errors = {};

%% Example 1: Search an events file for HED
request1 = struct('service', 'events_search', ...
                  'schema_version', '8.0.0', ...
                  'json_string', data.jsonText, ...
                  'events_string', data.eventsText, ...
                  'query', '[[Intended-effect, Cue]]');

response1 = webwrite(servicesUrl, request1, options);
response1 = jsondecode(response1);
outputReport(response1, 'Example 1 Querying an events file');
if ~isempty(response1.error_type) || ...
   ~strcmpi(response1.results.msg_category, 'success')
   errors{end + 1} = 'Example 1 failed execute the search.';
end

%% Example 2: Search an events file for HED
request2 = struct('service', 'events_search', ...
                  'schema_version', '8.0.0', ...
                  'json_string', data.jsonText, ...
                  'events_string', data.eventsText, ...
                  'columns_included', '', ...
                  'query', '[[Intended-effect, Cue]]');
request2.columns_included = {'onset'};
response2 = webwrite(servicesUrl, request2, options);
response2 = jsondecode(response2);
outputReport(response2, 'Example 2 Querying an events file with extra columns');
if ~isempty(response2.error_type) || ...
   ~strcmpi(response2.results.msg_category, 'success')
   errors{end + 1} = 'Example 2 failed execute the search.';
end

