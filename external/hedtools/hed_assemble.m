% Assemble HED tags for each event row given a tabular event file and a
% json sidecar containing HED annotations for event file columns.
%
% Usage:
%
%   >>  event = hed_assemble(events_tsv, events_json)
%
%
% Input:
%
%   Required:
%
%   events_tsv
%                    A tabular file containing event structure compatible 
%                    with BIDS events.tsv file
%
%   events_json
%                    A json sidecar containing HED annotations for columns
%                    in events_tsv in a format compatible with the BIDS
%                    events.json file
%
%
% Output:
%
%   evt_tbl
%                   An event struct array containing same columns and rows
%                   as in the events_tsv file with a column added
%                   containing assembled HED annotations for each event
%                   row. If an issue occurred, it returns an empty array.
%
% Copyright (C) 2022 Dung Truong dutruong@ucsd.edu,
% Kay Robbins kay.robbins@utsa.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
function evt_tbl = hed_assemble(events_tsv, events_json) 
    evt_tbl = [];
    jsonText = fileread(events_json);
    eventsText = fileread(events_tsv);
    host = 'https://hedtools.ucsd.edu/hed';
    csrfUrl = [host '/services']; 
    servicesUrl = [host '/services_submit'];
    [cookie, csrftoken] = getSessionInfo(csrfUrl);
    header = ["Content-Type" "application/json"; ...
              "Accept" "application/json"; ...
              "X-CSRFToken" csrftoken; "Cookie" cookie];

    options = weboptions('MediaType', 'application/json', 'Timeout', 120, ...
                         'HeaderFields', header);

    request = struct('service', 'events_assemble', ...
              'schema_version', '8.1.0', ...
              'json_string', jsonText, ...
              'events_string', eventsText, ...
              'expand_defs', 'on');
    
    try
        response = webwrite(servicesUrl, request, options);
    catch ME
        error('Issue connecting to HED tool web server. This could be due to internet connectivity. Please try again\n');
    end
    response = jsondecode(response);
    if ~isempty(response.error_type)
        error('Error assembling HED tags for the event file:\n\t-%s', response.error_msg);
    end

    if isfield(response, 'results') && ~isempty(response.results)
        results = response.results;
        fprintf('[%s] status %s: %s\n', response.service, results.msg_category, results.msg);
        fprintf('HED schema version: %s\n', results.schema_version);

        data = results.data;
        if ischar(data)
            data = strsplit(data, '\n');
        end
        if size(data,1) == 1 && size(data,2) > 1
            data = data';
        end
        HED = data(2:end-1); % first cell is header, last cell is empty
        hed_tbl = cell2table(HED);
        evt_tbl = readtable(events_tsv, 'FileType', 'delimitedtext', 'Delimiter', 'tab');
        if size(hed_tbl,1) == size(evt_tbl,1)
            evt_tbl = [evt_tbl hed_tbl];
        else
            error('Unexpected error: Assembled event table has different number of rows from that of original. Please report this issue to developer.');
        end
    else
        error('Unexpected error: empty results for successful assembling. Please report this issue to developer.');
    end
end