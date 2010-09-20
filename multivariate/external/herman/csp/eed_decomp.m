
function [F,d] = eed_decomp(filtered_data1,filtered_data2)
%
% EED_DECOMP is an auxiliary function used by csp_train (on its own or as
% part of CSPPROCESSOR object). It extracts spatial filtering matrix based
% on extreme energy difference criterion (EED) (see REFERENCES).
%
% Use as
%         1) Used by CSP_TRAIN (csp_train.m)
%         2) As part of CSPPROCESSOR object (see cspprocessor.m) - RECOMMENDED
%         3) [F,d] = eed_decomp(filtered_data1,filtered_data2)
%
% INPUT
%           filtered_data1, 
%           filtered_data2 - three-dimensional arrays [trial x chan x time] 
%                            of data corresponding to two classes
%
% REMARKS
%          1) Data can be prefiltered in the desirable frequency band (see PREPROCESSING)
%          2) Covariance is calculated along the 3rd dimension and averaged over the 1st one
%
% OUTPUT 
%           F  - arrays of EED filters derived from class1 and class2
%           d  - corresponding eigenvalues
%
% REFERENCES
%           Li J.,Sun S., Energy feature extraction of EEG signals and a case study, IJCNN 2008
%
% SEE ALSO
%           csp_test.m
%           csp_train.m
%           cspprocessor.m
%

% Pawel Herman, 2009

if size(filtered_data1,2) ~= size(filtered_data2,2)
    error('The second dim of the input arrays does not match');
else
    second_dim = size(filtered_data1,2); %the dimension for the cov calculations (often corresponding to channels)
end

data1.trial = num2cell(filtered_data1,[2,3])'; 
for i=1:length(data1.trial), data1.trial{i} = squeeze(data1.trial{i}); end
data1.dimord = 'rpt_chan_time'; 
data1.label = mat2cell(sprintf('%4d',1:second_dim),1,4*ones(1,second_dim))'; % it is not important what labels are attached to the second dim

data2.trial = num2cell(filtered_data2,[2,3])';
for i=1:length(data2.trial), data2.trial{i} = squeeze(data2.trial{i}); end
data2.dimord = 'rpt_chan_time'; 
data2.label = data1.label;

tlen = size(data1.trial{1},2);      
data1.time = mat2cell(repmat(1:tlen,1,length(data1.trial)),[1],tlen * ones(1,length(data1.trial)));
data2.time = mat2cell(repmat(1:tlen,1,length(data2.trial)),[1],tlen * ones(1,length(data2.trial)));

cfg.channel = data1.label;
cfg.keeptrials = 'no';         
cfg.covariance = 'yes';
cfg.covariancewindow  = [1 tlen];
cfg.normalizevar  = 'N-1';      
cfg.removemean = 'yes';
tlock1_struct = timelockanalysis(cfg,data1);
tlock2_struct = timelockanalysis(cfg,data2);
cov1 = tlock1_struct.cov;
cov2 = tlock2_struct.cov;

[F,d] = eig(cov1 - cov2);


