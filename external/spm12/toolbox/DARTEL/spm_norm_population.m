function out = spm_norm_population(job)
% Obtain mapping from population average to ICBM space
% FORMAT spm_norm_population(job)
% job.template - name of population average template
%
%________________________________________________________
% (c) Wellcome Trust Centre for NeuroImaging (2011)

% John Ashburner
% $Id: spm_norm_population.m 5677 2013-10-10 18:59:11Z john $

% Hard coded stuff, that should maybe be customisable
Ng  = nifti(fullfile(spm('Dir'),'toolbox','DARTEL','icbm152.nii'));
prm = [0 0.05 0.05 0.00000001 0.05 2 2 6 0];

if ~isempty(job.template{1})
    Nf = nifti(job.template{1});
else
    error('Incorrect Usage');
end

M0 = spm_klaff(Nf, Ng);    % Affine registration

df = [size(Nf.dat),1,1];   % Dimensions of data
dg = [size(Ng.dat),1,1];   % Dimensions of data

nd = min([df(4)+1,dg(4)]); % Use the first nd volumes.
nd = min(nd,2);            % Just use GM, WM & other

f  = single(Nf.dat(:,:,:,1:nd));   % Population average
g  = zeros([df(1:3),nd],'single'); % Resliced ICBM

% Images need to be same size, so reslice the MNI data
M = inv(M0);
for k=1:nd,
    p          = Ng.dat(:,:,:,k);
    [i1,i2,i3] = ndgrid(1:df(1),1:df(2),1:df(3));
    j1 = M(1,1)*i1+M(1,2)*i2+M(1,3)*i3+M(1,4);
    j2 = M(2,1)*i1+M(2,2)*i2+M(2,3)*i3+M(2,4);
    j3 = M(3,1)*i1+M(3,2)*i2+M(3,3)*i3+M(3,4);
    p  = spm_bsplins(p,j1,j2,j3,[1 1 1  0 0 0]);
    p(~isfinite(p)) = 0;
    g(:,:,:,k)      = p;
end

% Guess some settings that might work.  Note that this uses mean squared difference.
u   = zeros([df(1:3),3],'single'); % Starting estimates

spm_plot_convergence('Init','Least-squares Nonlin. Registration',...
              'Objective Fun.', 'Iteration');

prm1 = [prm(1) prm(2:4)/mean(sqrt(sum(Nf.mat(1:3,1:3).^2))) prm(5:end)];
for it=1:12,
    [u,ll] = dartel3(u,f,g,prm1); % Gauss-Newton update for registration
    fprintf('%d \t%g\t%g\t%g\t%g\n',...
        it,ll(1),ll(2),ll(1)+ll(2),ll(3));
    spm_plot_convergence('Set',ll(1)+ll(2));
end
spm_plot_convergence('Clear');

y1 = dartel3('Exp', u, [6,1]); % Mapping from voxels in f to voxels in g
M  = Ng.mat*inv(M0);           % Mapping from voxels in g to mm space of ICBM

% Compose the transforms
y2(:,:,:,1) = M(1,1)*y1(:,:,:,1)+M(1,2)*y1(:,:,:,2)+M(1,3)*y1(:,:,:,3)+M(1,4);
y2(:,:,:,2) = M(2,1)*y1(:,:,:,1)+M(2,2)*y1(:,:,:,2)+M(2,3)*y1(:,:,:,3)+M(2,4);
y2(:,:,:,3) = M(3,1)*y1(:,:,:,1)+M(3,2)*y1(:,:,:,2)+M(3,3)*y1(:,:,:,3)+M(3,4);

% Set up output deformation header
[pth,nam,ext]       = fileparts(job.template{1});
defname             = fullfile(pth,['y_' nam '_2mni.nii']);
dat                 = file_array(defname,[df(1:3),1,3],'float32',0,1,0);
Nout                = nifti;
Nout.dat            = dat;
Nout.mat            = Nf.mat;
Nout.mat_intent     = 2;
Nout.descrip        = 'Population to MNI';

% Write output deformation
create(Nout);
Nout.dat(:,:,:,1,1) = y2(:,:,:,1);
Nout.dat(:,:,:,1,2) = y2(:,:,:,2);
Nout.dat(:,:,:,1,3) = y2(:,:,:,3);

out.files{1} = defname;

