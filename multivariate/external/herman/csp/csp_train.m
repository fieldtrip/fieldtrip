
function [filters,d] = csp_train(data,design,k_patterns,mode)
%
% CSP_TRAIN extracts CSP/EED filters from two-class data
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
%                         TIMELOCKANALYSIS(append([],preproc1,preproc2)) 
%                         (design vector has to be created and appropriately concatenated to match assigments of data to conditions)
%                   2) 3-D array of rpt_chan_time format - e.g. as a result of concatenation of trial arrays obtained using 
%                      TIMELOCKANALYSIS -> [tlock1.trial; tlock2.trial] 
%                         (design vector has to be created and appropriately concatenated to match assigments of data to conditions)
%
%           design     - design matrix referring to class labels associated
%                        with input data (dichotomous problems)
%
%           k_patterns - number of CSP/EED filters/patterns to be extracted
%
%           mode       - type of algorithm to be used for spatial filtering
%                         'CSP0' -> CSP by common diagonalisation of class1+class2 
%                                  covariance matrix (biosig-> csp.m, see also csp_wrapper.m)
%                         'CSP3' -> CSP calculation as generalized eigenvalues
%                                   (biosig-> csp.m, see also csp_wrapper.m)
%                       'stdCOV' -> based on standard Fieldtrip estimation of
%                                   covariance matrix 
%                          'EED' -> extreme energy difference criterion
%
%
% REMARKS
%           All three input data formats can be used when csp_train is
%           called outside typical CLFPROC pipeline. Otherwise, only 3D
%           array representation is allowed.
%
% OUTPUT 
%           filters - arrays of CSP filters (2 x k_patterns column vectors)
%           d       - corresponding eigenvalues (in pairs - the first one is the
%                     largest for class1 and the last the smallest for class2)
%
% SEE ALSO
%           csp_test.m
%           csp_train.m
%           cspprocessor.m
%

% Pawel Herman, 2009

if nargin<3, k_patterns = 1; end
if nargin<4, mode = 'fieldtrip_cov'; end

labels = unique(design);
if length(labels)~=2,  error('For the simplified CSP calculation, two-class data set is needed.'); end

if isstruct(data)
    if isfield(data,'trial')
        if iscell(data.trial)
            ind1 = find(design==labels(1));
            ind2 = find(design==labels(2));            
            filt_data1 = reshape([data.trial{ind1}],[size(data.trial{ind1(1)}),length(ind1)]);
            filt_data1 = permute(filt_data1,[3 1 2]);
            filt_data2 = reshape([data.trial{ind2}],[size(data.trial{ind2(1)}),length(ind2)]);
            filt_data2 = permute(filt_data2,[3 1 2]);            
        elseif isnumeric(data.trial)
            filt_data1 = data.trial(design==labels(1),:,:);
            filt_data2 = data.trial(design==labels(2),:,:);            
        else
            error('Input data is a structure with unrecognized format of field .trial');
        end
    else
        error('Input data is a structure without field .trial');
    end
elseif isnumeric(data)
    filt_data1 = data(design==labels(1),:,:);
    filt_data2 = data(design==labels(2),:,:);
else
    error('Unrecognized data format');
end

if strcmp(mode,'biosig_CSP0') || strcmp(mode,'CSP0')
    [filters,d] = csp_wrapper(filt_data1,filt_data2,k_patterns,'CSP0');
    filters = filters(:,1:end/2);
    d = d(1:end/2);    
    
elseif strcmp(mode,'biosig_CSP3') || strcmp(mode,'CSP3')
    [filters,d] = csp_wrapper(filt_data1,filt_data2,k_patterns,'CSP3');
    filters = filters(:,1:end/2);
    d = d(1:end/2);
    
elseif strcmp(mode,'stdCOV')
    [V1,d1,V2,d2] = csp_decomp(filt_data1,filt_data2);
    len1 = length(d1); len2 = length(d2);   
    assert(len1==len2);
    
    filters = V1([1:k_patterns len1-k_patterns+1:len1],:)';
    %filters = [filters V2([1:k_patterns len2-k_patterns+1:len2],:)'];
    d = d1([1:k_patterns len1-k_patterns+1:len1])';
    %d = [ d d2([1:k_patterns len2-k_patterns+1:len2])' ];   
    
elseif strcmp(mode,'EED')  %extreme energy difference criterion
    [filters,d] = eed_decomp(filt_data1,filt_data2);
    len = length(d);
    filters = filters([1:k_patterns len-k_patterns+1:len],:)'; 
    d = d([1:k_patterns len-k_patterns+1:len])';
end


