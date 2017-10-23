function [q,res] = spm_mars_core(q,res,cr,tpm,J,lkp,rLabel,convergence)
% [q,res] = spm_mars_core(q,res,cr,tpm,J,lkp,rLabel,convergence)
%
% Core function of the toolbox
% 'MARS (Morphologically and Anatomically accuRate Segmentation)'.
%
% EM iterations are integrated with MRF constraints. Checkerboard updating scheme.
%
% q: input posterior probabilities, and also output posterior probabilities after iteration
% res: a structure containing all related info from the results of 'New Segment'
% cr: bias corrected image(s)
% tpm: tissue probability map (mean Mi = <xi> for each voxel)
% J: interaction term from TCM (global, local or regional)
% lkp: for compatability with 'New Segment', i.e., the multi-gaussian for one
% tissue type; also compatible with non-parametric intensity model
% rLabel: regional labels, for regional TCM only
% convergence: user-defined threshold for the convergence of the algorithm
%
% Some codes are adapted from spm_preproc8.m by John Ashburner
%
% Yu (Andy) Huang, 2013-05-15
% Yu (Andy) Huang, 2014-08-30
% $Id: spm_mars_core.m 2015-07-27 andy$
% Neural Engineering Lab, Dept. of Biomedical Engineering, City College of New York
% yhuang16@citymail.cuny.edu

tiny = 1.0e-64; % andy 2013-07-24 % USE SUPER SMALL VALUE. IF THIS IS NOT SMALL ENOUGH, THEN IT WILL CAUSE NUMERICAL ERROR. BUT IT CANNOT BE TOO SMALL.

% Compute the total F of the image %andy 2012-12-20
Kb = max(lkp);
K = numel(lkp);
N = numel(res.image);
d = res.image(1).dim(1:3);
nm = prod(d);

%================convergence study  % andy 2013-07-03======================
oV = zeros(Kb,1);
for k=1:Kb
    tmp = sum(q(:,:,:,lkp==k),4);
    oV(k) = sum(tmp(:));
end
eV = 1;
%================convergence study  % andy 2013-07-03======================

ext_eng = log(tpm+tiny); % external energy from TPM  % andy 2013-03-27
lkp0 = int32(lkp-1); % just for C routines, since index starts from 0 in C, and 'int' in C has 32bits  % andy 2013-05-15

cr = reshape(cr,[nm N]);

FEF = 0; % for fast speed, it is NOT necessary to compute an initial total free energy, because ONLY the change
% of FEF is needed to terminate the EM iterations, so set initial FEF to 0. % andy 2013-05-15

if isfield(res,'mg')
    % Add a little something to the covariance estimates in order to assure stability
    vr0 = zeros(N,N);
    for n=1:N,
        %         if spm_type(V(n).dt(1),'intt'),
        vr0(n,n) = 0.083; %*V(n).pinfo(1,1);
        %         else
        %             vr0(n,n) = (max(mn(n,:))-min(mn(n,:)))^2*1e-8;
        %         end
    end % NEED TO BE ADDED IN THE FINAL VERSION, NOT ALWAYS USE 0.083  % andy 2013-04-30
else
    intensity_tmp = zeros(nm,N);
    for n=1:N
        tmp = round(cr(:,n)*res.intensity(n).interscal(2) + res.intensity(n).interscal(1));
        intensity_tmp(:,n) = min(max(tmp,1),size(res.intensity(n).lik,1));
        % modeled (not observed) image intensities used for estimating
        % histogram for each tissue type
    end
    x = (1:size(res.intensity(1).lik,1))';
end

spm_plot_convergence('Init','Processing','Total Free Energy (Initial Value Set to 0)','Iteration');
%------------------------------------------------------------
subit = 1;
% fid = fopen('convergence','w'); % andy 2013-05-04
% while eV >= 1e-04 % use volume change as convergence criterion % andy 2013-08-02
while eV >= convergence
%     ts = tic;
    fprintf('Iteration #%d...\n',subit);
    %     fprintf(fid,'Iteration #%d...\r\n',subit); % andy 2013-05-04
%     fprintf(fid,'Iteration #%d...\n',subit); % andy 2013-06-12
    
    if isfield(res,'mg')
        l = likelihoods(cr,res.mg,res.mn,res.vr);
        l = reshape(l,[d(1:3) K]); % likelihoods computed for all voxels at one time % andy 2013-03-27
    else
        l = ones(nm,Kb);
        for n=1:N
            for k=1:Kb,
                likelihood = res.intensity(n).lik(:,k);
                l(:,k) = l(:,k) .* likelihood(intensity_tmp(:,n));
            end
        end
        l = reshape(l,[d(1:3) Kb]);
    end
    
    %-----------------------------------------------------------------------
    % Asynchronous updating for posterior q_ik
    % CORE UPDATING STEP (E-STEP) FOR THE POSTERIOR q_ik, C routine % andy 2013-02-22
    deltaFEF = spm_mars(J,q,l,ext_eng,lkp0,rLabel);
    %-----------------------------------------------------------------E STEP
    
    if isfield(res,'mg')
        mom0 = zeros(K,1)+tiny;
        mom1 = zeros(N,K);
        mom2 = zeros(N,N,K);
        
        for k=1:K, % Moments
            tmp = reshape(q(:,:,:,k),nm,1);
            mom0(k)     = sum(tmp);
            mom1(:,k)   = (tmp'*cr)'; % (tmp'*reshape(cr,nm,N))';
            mom2(:,:,k) = (repmat(tmp,1,N).*cr)'*cr; % (repmat(tmp,1,N).*reshape(cr,nm,N))'*reshape(cr,nm,N);
        end
        
        % Mixing proportions, Means and Variances
        for k=1:K,
            tmp       = mom0(lkp==lkp(k));
            res.mg(k)     = (mom0(k)+tiny)/sum(tmp+tiny);
            res.mn(:,k)   = mom1(:,k)/(mom0(k)+tiny);
            res.vr(:,:,k) = (mom2(:,:,k) - mom1(:,k)*mom1(:,k)'/mom0(k))/(mom0(k)+tiny) + vr0; %vr0 for stability
        end % update equations for mu, sigma, and gamma
        % NOTE: mu(mn), sigma(vr) are both consistent with the math formula in the paper; ONLY gamma(mg) is different from paper (mentioned at beginning of spm_preproc8). andy 2013-04-30
        
    else
        for n=1:N,
            h = zeros(size(res.intensity(n).lik,1),Kb);
            lam = zeros(Kb,1);
            for k=1:Kb,
                h(:,k) = accumarray(intensity_tmp(:,n),reshape(q(:,:,:,k),nm,1),[size(res.intensity(n).lik,1),1]);
                % historgram
            end
            
            for k=1:Kb,
                mom0 = sum(h(:,k)) + eps;
                mom1 = sum(h(:,k).*x) + eps;
                lam(k) = sum(h(:,k).*(x-mom1./mom0).^2+1)/(mom0+1)+1;
            end
            
            res.intensity(n).lik   = spm_smohist(h,lam)*res.intensity(n).interscal(2);
            % from histogram to pdf
        end
        
    end
    %-----------------------------------------------------------------M STEP
    
    FEF = FEF + deltaFEF; % andy 2013-05-14
    spm_plot_convergence('Set',FEF); %andy 2012-12-20
    
%     fprintf('Change of total free energy after updating all voxels is %f.\n',deltaFEF); % andy 2013-05-14
    %     fprintf('Change of total free energy after updating all voxels is %f.\n',FEF - oFEF); %andy 2012-12-20
    %     fprintf(fid,'Change of total free energy after updating all voxels is %f.\r\n',FEF - oFEF); % andy 2013-05-04
%     fprintf(fid,'Change of total free energy after updating all voxels is %f.\n',deltaFEF); % andy 2013-06-12
    
%     toc(ts);
    
    %================convergence study  % andy 2013-07-03======================
    V = zeros(Kb,1);
    for k=1:Kb
        tmp = sum(q(:,:,:,lkp==k),4);
        V(k) = sum(tmp(:));
    end % convergence study  % andy 2013-07-03
    deltaV = V - oV;
%     fprintf('Change of posterior volume after updating all voxels is:\n%f\t%f\t%f\t%f\t%f\t%f\n',deltaV(1),deltaV(2),deltaV(3),deltaV(4),deltaV(5),deltaV(6));
%     fprintf(fid,'Change of posterior volume after updating all voxels is:\n%f\t%f\t%f\t%f\t%f\t%f\n',deltaV(1),deltaV(2),deltaV(3),deltaV(4),deltaV(5),deltaV(6));
    rdeltaV = deltaV./oV;
%     fprintf('Relative change of posterior volume is:\n%f\t%f\t%f\t%f\t%f\t%f\n',rdeltaV(1),rdeltaV(2),rdeltaV(3),rdeltaV(4),rdeltaV(5),rdeltaV(6));
%     fprintf(fid,'Relative change of posterior volume is:\n%f\t%f\t%f\t%f\t%f\t%f\n',rdeltaV(1),rdeltaV(2),rdeltaV(3),rdeltaV(4),rdeltaV(5),rdeltaV(6));
    eV = max(abs(rdeltaV));
    oV = V;
    %================convergence study  % andy 2013-07-03======================
    subit = subit + 1;
    
    if subit>=60, break; end % iterations more than 60 are meaningless % andy 2014-03-31
end
% fclose(fid); % andy 2013-05-04

fprintf('MARS iterations done!\n');
res.FEF      = FEF;

function l = likelihoods(intensity,mg,mn,vr)
K  = numel(mg);
% N  = numel(intensity);
[M N] = size(intensity);
l  = zeros(M,K);
for k=1:K,
    amp    = mg(k)/sqrt((2*pi)^N * det(vr(:,:,k)));
    d      = intensity - repmat(mn(:,k)',M,1);
    l(:,k) = amp * exp(-0.5* sum(d.*(d/vr(:,:,k)),2));
    % -0.5* sum(d(i,:).*(d(i,:)/vr(:,:,k)),2) is equivalent to -0.5 * d(i,:) * inv(vr(:,:,k)) * d(i,:)'
    % i stands for one location
end
% make likelihoods function be able to compute likelihoods for the entire
% image % andy 2013-03-25