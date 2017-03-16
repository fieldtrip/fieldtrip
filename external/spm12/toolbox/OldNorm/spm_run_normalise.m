function out = spm_run_normalise(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_normalise.m 4873 2012-08-30 19:06:26Z john $


if isfield(job,'eoptions')
    eflags = rmfield(job.eoptions,{'weight','template'});
end
    
for i=1:numel(job.subj)
    
    %-Deformation file
    %----------------------------------------------------------------------
    if isfield(job.subj(i),'matname')
        params = char(job.subj(i).matname);
    else
        params = spm_file(char(job.subj(i).source), 'suffix','_sn', 'ext','.mat');
    end
        
    %-Normalise: Estimate
    %----------------------------------------------------------------------
    if isfield(job,'eoptions')
        
        spm_normalise(char(job.eoptions.template),...
            char(job.subj(i).source), params,...
            char(job.eoptions.weight), char(job.subj(i).wtsrc), eflags);
    end
    
    %-Normalise: Write
    %----------------------------------------------------------------------
    if isfield(job,'roptions')
        spm_write_sn(char(job.subj(i).resample), params, job.roptions);
    end
end

%-Dependencies
%--------------------------------------------------------------------------
for i=1:numel(job.subj)
    if ~isfield(job.subj(i),'matname')
        out(i).params = {spm_file(char(job.subj(i).source), 'suffix','_sn', 'ext','.mat')};
    end
    
    if isfield(job,'roptions')
        out(i).files = spm_file(job.subj(i).resample, 'prefix',job.roptions.prefix);
    end
end
