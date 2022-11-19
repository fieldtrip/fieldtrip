function errors = testGetServices(host)
%% Get the options and data
[servicesUrl, options] = getHostOptions(host);
errors = {};

%% Send the request and get the response
request = struct('service', 'get_services');
response = webwrite(servicesUrl, request, options);
response = jsondecode(response);
fprintf('Error report:  [%s] %s\n', response.error_type, response.error_msg);

%% Print out the results if available
if isfield(response, 'results') && ~isempty(response.results)
   results = response.results;
   fprintf('[%s] status %s: %s\n', response.service, results.msg_category, results.msg);
   fprintf('Return data:\n%s\n', results.data);
end

if ~isempty(response.error_type) || ...
   ~strcmpi(response.results.msg_category, 'success')
   errors{end + 1} = 'Get services failed to return services.';
end