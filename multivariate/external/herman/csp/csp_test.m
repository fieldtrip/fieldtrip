
function [csp_proj,csp_pow] = csp_test(data,filters)

%
% CSP_TEST applies CSP/EED filters to data
%
% Use as
%         1) [filters,d] = csp_train(data,design,k_patterns,mode)
%         2) As part of CSPPROCESSOR object (see cspprocessor.m)
%
% INPUT
%           data  - input data can have three different formats (1a, 1b or 2):                       
%                   1) structure with field .trial 
%                      a) trial can be a cell array with chan_time data representation 
%                         (the result of timelockanalysis and append([],tlock1,tlock2) 
%                         for two conditions)  OR 
%                      b) trial can be a 3-D array of rpt_chan_time format - 
%                         concatenation of trial arrays obtained using
%                         TIMELOCKANALYSIS -> [tlock1.trial; tlock2.trial] or as a result of 
%                         TIMELOCKANALYSIS(append([],preproc1,preproc2)) %                         
%                   2) 3-D array of rpt_chan_time format - e.g. as a result of concatenation of trial arrays obtained using 
%                      TIMELOCKANALYSIS -> [tlock1.trial; tlock2.trial] 
%
%           filters - array of filters (represented as column vectors)
%
%
% REMARKS
%           All three input data formats can be used when csp_test is
%           called outside typical CLFPROC pipeline. Otherwise, only 3D
%           array representation is allowed.
%
% OUTPUT 
%           csp_proj  - spatially filtered signals (CSP or EED projection)
%           csp_pow   - power of the projected (spatially filtered) signals
%
% SEE ALSO
%           csp_test.m
%           csp_train.m
%           cspprocessor.m
%

% Pawel Herman, 2009

csp_pow = [];

if isstruct(data)
    if isfield(data,'trial')
        if iscell(data.trial)
            filt_data = permute(reshape([data.trial{:}],[size(data.trial{1}),length(data.trial)]),[3 1 2]);
            dat_option = '1a';  
            csp_proj = data; csp_proj.trial = {};
        elseif isnumeric(data.trial)
            filt_data = data.trial;
            dat_option = '1b'; 
            csp_proj = data; csp_proj.trial = [];
        else
            error('Input data is a structure with unrecognized format of field .trial');
        end
    else
        error('Input data is a structure without field .trial');
    end
elseif isnumeric(data)
    filt_data = data;
    dat_option = '2';   csp_proj = [];
else
    error('Unrecognized data format');
end

for i=1:size(filt_data,1)
    switch dat_option
        case '1a',
            csp_proj.trial{i} = filters' * squeeze(filt_data(i,:,:));
            csp_pow = [csp_pow;  (sum(csp_proj.trial{i}.^2,2)')];
        case '1b',
            csp_proj.trial(i,:,:) = filters' * squeeze(filt_data(i,:,:));
            csp_pow = [csp_pow;  (sum(squeeze(csp_proj.trial(i,:,:)).^2,2)')];
        case '2',
            csp_proj(i,:,:) = filters' * squeeze(filt_data(i,:,:));
            csp_pow = [csp_pow;  (sum(squeeze(csp_proj(i,:,:)).^2,2)')];
    end

end



