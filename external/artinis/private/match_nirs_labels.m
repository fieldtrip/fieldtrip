function matching_channels = match_nirs_labels(labels)
% MATCHING_CHANNELS finds the labels of the signals associated to the
% same combination of receiver and transmitter and provides as output a
% binary matrix where 1 represents that channels in the corresponding
% row-column combination match.
% 
% OUTPUT
%  matching_channels: nChan x nChan binary matrix; where 1 represents that 
% channels in the corresponding row-column combination match, and 0 represents
% the cases where they do not match.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_channels=numel(labels);
matching_channels = ones(num_channels, num_channels);

for n = 1 : num_channels
    label_check1 = labels{n};
    label_check1_split = split(label_check1, ' [');
    label_check1_split = label_check1_split{1};
    for m = 1 : num_channels
        if(n == m)
            continue
        end
        label_check2 = labels{m};
        label_check2_split = split(label_check2, ' [');
        label_check2_split = label_check2_split{1};
        matching_channels(n,m) = strcmp(label_check1_split, label_check2_split);
    end
end



    
    
    
    
    
    