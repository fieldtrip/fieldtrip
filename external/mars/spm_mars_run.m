function varargout = spm_mars_run(job,arg)
% Entry function of the toolbox
% 'MARS (Morphologically and Anatomically accuRate Segmentation)'
%
% FORMAT spm_mars_run(job)
%
% job.channel(n).vols{m}  (n is index for channels, e.g. T1, T2, etc.; m is
% index for subjects)
% job.channel(n).biasreg
% job.channel(n).biasfwhm
% job.channel(n).write
% job.tissue(k).tpm  (k is index for tissue types)
% job.tissue(k).ngaus
% job.tissue(k).native
% job.tissue(k).warped
% job.warp.reg
% job.warp.affreg
% job.warp.samp
% job.warp.write
% job.marsOptions.marsMode
% job.marsOptions.tcm
% job.marsOptions.beta
% job.marsOptions.marsConvergence
%
% See the user interface for a description of the fields.
%
% Adapted from spm_preproc_run (entry function for toolbox 'New Segment') by
% John Ashburner
% $Id: spm_preproc_run.m 4334 2011-05-31 16:39:53Z john $
% Copyright (C) 2008-2011 Wellcome Trust Centre for Neuroimaging
%
% Yu (Andy) Huang
% $Id: spm_mars_run.m 2015-07-16 andy$
% Neural Engineering Lab, Dept. of Biomedical Engineering, City College of New York
% yhuang16@citymail.cuny.edu

% if nargin == 1, arg = 'run'; end

% switch lower(arg)
%     case 'run'
        varargout{1} = run_job(job);
%     case 'check'
%         varargout{1} = check_job(job);
%     case 'vfiles'
%         varargout{1} = vfiles_job(job);
%     case 'vout'
%         varargout{1} = vout_job(job);
%     otherwise
%         error('Unknown argument ("%s").', arg);
% end
return
%_______________________________________________________________________

%_______________________________________________________________________
function vout = run_job(job)

vout   = vout_job(job);
tpm    = strvcat(cat(1,job.tissue(:).tpm));
tpm    = spm_load_priors8(tpm);

% nit = 1;

% for iter=1:nit,
%     if nit>1,
%         % Sufficient statistics for possible generation of group-specific
%         % template data.
%         SS = zeros([size(tpm.dat{1}),numel(tpm.dat)],'single');
%     end
for subj=1:numel(job.channel(1).vols), % for multiple subjects
    images = '';
    for n=1:numel(job.channel),
        images = strvcat(images,job.channel(n).vols{subj});
    end
    obj.image    = spm_vol(images);
    spm_check_orientations(obj.image);
    
    obj.fudge    = 5;
    obj.biasreg  = cat(1,job.channel(:).biasreg);
    obj.biasfwhm = cat(1,job.channel(:).biasfwhm);
    obj.tpm      = tpm;
    obj.lkp      = [];
    if all(isfinite(cat(1,job.tissue.ngaus))),
        for k=1:numel(job.tissue),
            obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
        end;
    end
    obj.reg      = job.warp.reg;
    obj.samp     = job.warp.samp;
    
    %             if iter==1,
    % Initial affine registration.
    Affine  = eye(4);
    if ~isempty(job.warp.affreg),
        Affine  = spm_maff8(obj.image(1),job.warp.samp,obj.fudge*8,tpm,Affine,job.warp.affreg); % Close to rigid
        Affine  = spm_maff8(obj.image(1),job.warp.samp,obj.fudge,  tpm,Affine,job.warp.affreg);
    end;
    obj.Affine = Affine;
    %             else
    %                 % Load results from previous iteration for use with next round of
    %                 % iterations, with the new group-specific tissue probability map.
    %                 [pth,nam] = fileparts(job.channel(1).vols{subj});
    %                 res       = load(fullfile(pth,[nam '_seg8.mat']));
    %                 obj.Affine = res.Affine;
    %                 obj.Twarp  = res.Twarp;
    %                 obj.Tbias  = res.Tbias;
    %                 if ~isempty(obj.lkp),
    %                     obj.mg     = res.mg;
    %                     obj.mn     = res.mn;
    %                     obj.vr     = res.vr;
    %                 end
    %             end
    
    res = spm_preproc8(obj);
    
    try
        [pth,nam] = fileparts(job.channel(1).vols{subj});
        savefields(fullfile(pth,[nam '_seg8.mat']),res);
    catch
    end
    
    %         if iter==nit,
    % Final iteration, so write out the required data.
    tmp1 = [cat(1,job.tissue(:).native) cat(1,job.tissue(:).warped)];
    tmp2 =  cat(1,job.channel(:).write);
    tmp3 = job.warp.write;
    if job.marsOptions.marsMode == 0
        fprintf('New Segment done, writing out results without post-processing...\n');
        spm_mars_newSeg(res,tmp1,tmp2,tmp3,job.marsOptions.marsMode,job.marsOptions.tcm);
    elseif job.marsOptions.marsMode == 1
        fprintf('New Segment done, writing out results with MRF cleanup as post-processing...\n');
        spm_mars_newSeg(res,tmp1,tmp2,tmp3,job.marsOptions.marsMode,job.marsOptions.tcm);
        % New Segment (without/with MRF cleanup) output results % andy 2013-05-03
    elseif job.marsOptions.marsMode == 2
        fprintf('New Segment done, begin to run MARS...\n');
        res = spm_mars_mrf(res,tmp1,tmp2,tmp3,job.marsOptions.tcm,job.marsOptions.beta,job.marsOptions.marsConvergence);
        % MARS initialized by New Segment % andy 2013-03-06
        
        try
            [pth,nam] = fileparts(job.channel(1).vols{subj}); % save latest parameters
            savefields(fullfile(pth,[nam '_seg8.mat']),res);
        catch
        end
    end
    %         else
    %             % Not the final iteration, so compute sufficient statistics for
    %             % re-estimating the template data.
    %             N    = numel(job.channel);
    %             K    = numel(job.tissue);
    %             cls  = spm_preproc_write8(res,zeros(K,4),zeros(N,2),[0 0],job.warp.mrf);
    %             for k=1:K,
    %                 SS(:,:,:,k) = SS(:,:,:,k) + cls{k};
    %             end
    %         end
end
%     if iter<nit && nit>1,
%         % Treat the tissue probability maps as Dirichlet priors, and compute the
%         % MAP estimate of group tissue probability map using the sufficient
%         % statistics.
%         alpha = 1.0;
%         for k=1:K,
%             SS(:,:,:,k) = SS(:,:,:,k) + spm_bsplinc(tpm.V(k),[0 0 0  0 0 0])*alpha + eps;
%         end
%         %save SS.mat SS
%         s = sum(SS,4);
%         for k=1:K,
%             tmp        = SS(:,:,:,k)./s;
%             tpm.bg1(k) = mean(mean(tmp(:,:,1)));
%             tpm.bg2(k) = mean(mean(tmp(:,:,end)));
%             tpm.dat{k} = spm_bsplinc(log(tmp+tpm.tiny),[ones(1,3)*(tpm.deg-1)  0 0 0]);
%         end
%     end
% end
return
%_______________________________________________________________________

% %_______________________________________________________________________
% function msg = check_job(job)
% msg = {};
% if numel(job.channel) >1,
%     k = numel(job.channel(1).vols);
%     for i=2:numel(job.channel),
%         if numel(job.channel(i).vols)~=k,
%             msg = {['Incompatible number of images in channel ' num2str(i)]};
%             break
%         end
%     end
% elseif numel(job.channel)==0,
%     msg = {'No data'};
% end
% return
% %_______________________________________________________________________
% 
%_______________________________________________________________________
function savefields(fnam,p)
if length(p)>1, error('Can''t save fields.'); end
fn = fieldnames(p);
if numel(fn)==0, return; end
for i=1:length(fn)
    eval([fn{i} '= p.' fn{i} ';']);
end
if spm_check_version('matlab','7') >= 0
    save(fnam,'-V6',fn{:});
else
    save(fnam,fn{:});
end

return;
%_______________________________________________________________________

%_______________________________________________________________________
function vout = vout_job(job)

n     = numel(job.channel(1).vols);
parts = cell(n,4);

channel = struct('biasfield',{},'biascorr',{});
for i=1:numel(job.channel),
    for j=1:n,
        [parts{j,:}] = spm_fileparts(job.channel(i).vols{j});
    end
    if job.channel(i).write(1),
        channel(i).biasfield = cell(n,1);
        for j=1:n
            channel(i).biasfield{j} = fullfile(parts{j,1},['BiasField_',parts{j,2},'.nii']);
        end
    end
    if job.channel(i).write(2),
        channel(i).biascorr = cell(n,1);
        for j=1:n
            channel(i).biascorr{j} = fullfile(parts{j,1},['m',parts{j,2},'.nii']);
        end
    end
end

for j=1:n,
    [parts{j,:}] = spm_fileparts(job.channel(1).vols{j});
end
param = cell(n,1);
for j=1:n
    param{j} = fullfile(parts{j,1},[parts{j,2},'_seg8.mat']);
end

tiss = struct('c',{},'rc',{},'wc',{},'mwc',{});
for i=1:numel(job.tissue),
    if job.tissue(i).native(1),
        tiss(i).c = cell(n,1);
        for j=1:n
            tiss(i).c{j} = fullfile(parts{j,1},['c',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).native(2),
        tiss(i).rc = cell(n,1);
        for j=1:n
            tiss(i).rc{j} = fullfile(parts{j,1},['rc',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).warped(1),
        tiss(i).wc = cell(n,1);
        for j=1:n
            tiss(i).wc{j} = fullfile(parts{j,1},['wc',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).warped(2),
        tiss(i).mwc = cell(n,1);
        for j=1:n
            tiss(i).mwc{j} = fullfile(parts{j,1},['mwc',num2str(i),parts{j,2},'.nii']);
        end
    end
end

if job.warp.write(1),
    invdef = cell(n,1);
    for j=1:n
        invdef{j} = fullfile(parts{j,1},['iy_',parts{j,2},'.nii']);
    end
else
    invdef = {};
end

if job.warp.write(2),
    fordef = cell(n,1);
    for j=1:n
        fordef{j} = fullfile(parts{j,1},['y_',parts{j,2},'.nii']);
    end
else
    fordef = {};
end

vout  = struct('channel',channel,'tiss',tiss,'param',{param},'invdef',{invdef},'fordef',{fordef});
%_______________________________________________________________________

% %_______________________________________________________________________
% function vf = vfiles_job(job)
% vout = vout_job(job);
% vf   = vout.param;
% if ~isempty(vout.invdef), vf = {vf{:}, vout.invdef{:}}; end
% if ~isempty(vout.fordef), vf = {vf{:}, vout.fordef{:}}; end
% for i=1:numel(vout.channel),
%     if ~isempty(vout.channel(i).biasfield), vf = {vf{:}, vout.channel(i).biasfield{:}}; end
%     if ~isempty(vout.channel(i).biascorr),  vf = {vf{:}, vout.channel(i).biascorr{:}};  end
% end
% 
% for i=1:numel(vout.tiss)
%     if ~isempty(vout.tiss(i).c),   vf = {vf{:}, vout.tiss(i).c{:}};   end
%     if ~isempty(vout.tiss(i).rc),  vf = {vf{:}, vout.tiss(i).rc{:}};  end
%     if ~isempty(vout.tiss(i).wc),  vf = {vf{:}, vout.tiss(i).wc{:}};  end
%     if ~isempty(vout.tiss(i).mwc), vf = {vf{:}, vout.tiss(i).mwc{:}}; end
% end
% vf = reshape(vf,numel(vf),1);
% %_______________________________________________________________________
% 
% %_______________________________________________________________________
