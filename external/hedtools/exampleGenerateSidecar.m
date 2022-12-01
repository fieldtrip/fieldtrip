host = 'https://hedtools.ucsd.edu/hed';
[servicesUrl, options] = getHostOptions(host); % Setup the options
pathname = '../../../datasets/eeg_ds003654s_hed/sub-002/eeg/';
filename = 'sub-002_task-FacePerception_run-1_events.tsv';
eventsText = fileread([pathname filename]); % Read the data
request = struct('service', 'events_generate_sidecar', ...
                 'events_string', eventsText);
request.columns_categorical = {'event_type', 'face_type', 'rep_status'};
request.columns_value = {'trial', 'rep_lag', 'stim_file'};
response = webwrite(servicesUrl, request, options);
response = jsondecode(response);
outputReport(response,'Example: generate a sidecar from an event file.');
