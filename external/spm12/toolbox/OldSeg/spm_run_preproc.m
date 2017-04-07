function out = spm_run_preproc(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_preproc.m 4873 2012-08-30 19:06:26Z john $

job.opts.tpm = char(job.opts.tpm);
job.opts.msk = char(job.opts.msk);

for i=1:numel(job.data)
    %-Combined Segmentation and Spatial Normalisation
    %----------------------------------------------------------------------
    res = spm_preproc(job.data{i}, job.opts);
    
    %-Convert the output from previous step into an sn.mat file
    %----------------------------------------------------------------------
    spm_prep2sn(res);
    snfile{i} = spm_file(job.data{i}, 'suffix','_seg_sn', 'ext','.mat');
end

%-Write out preprocessed data
%--------------------------------------------------------------------------
spm_preproc_write(snfile, job.output);

%-Dependencies
%--------------------------------------------------------------------------
for i=1:numel(job.data)
    out.snfile{i} = spm_file(job.data{i}, 'suffix','_seg_sn', 'ext','.mat');
    out.isnfile{i} = spm_file(job.data{i}, 'suffix','_seg_inv_sn', 'ext','.mat');
end

opts  = job.output;
sopts = [opts.GM;opts.WM;opts.CSF];
for i=1:numel(job.data)
    if opts.biascor
        out.biascorr{i,1} = spm_file(job.data{i}, 'prefix','m');
    end
    for k=1:3
        if sopts(k,3)
            out.(sprintf('c%d',k)){i,1} = ...
                spm_file(job.data{i}, 'prefix',sprintf('c%d',k));
        end
        if sopts(k,2)
            out.(sprintf('wc%d',k)){i,1} = ...
                spm_file(job.data{i}, 'prefix',sprintf('wc%d',k));
        end
        if sopts(k,1)
            out.(sprintf('mwc%d',k)){i,1} = ...
                spm_file(job.data{i}, 'prefix',sprintf('mwc%d',k));
        end
    end
end
