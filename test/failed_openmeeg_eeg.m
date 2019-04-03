function test_openmeeg_eeg

% WALLTIME 00:30:00
% MEM 4gb

% Test the computation of an EEG leadfield with OpenMEEG
% Copyright (C) 2010-2017, OpenMEEG developers

addpath(pwd) % Make sure current folder is in the path

%% Set the position of the probe dipole
dippos = [0 0 70];

%% Set the radius and conductivities of each of the compartments

% 4 Layers, defined from outside to inside (i.e., [scalp skull CSF brain])
r = [100 92 88 85];
c = [1 1/80 1/20 1];

[rdms,mags] = run_bem_computation(r,c,dippos);

% the following would require the installation of xunit toolbox
% assertElementsAlmostEqual(rdms, [0.019963 0.019962 0.10754], 'absolute', 1e-3)
% assertElementsAlmostEqual(mags, [0.84467 0.84469 0.83887], 'absolute', 1e-3)

%use instead
thr = 2e-2;
assert(norm(rdms-[0.019963 0.019962 0.10754])<thr)
assert(norm(mags-[0.84467 0.84469 0.83887])<thr)

% 3 Layers, defined from outside to inside (i.e., [scalp skull brain])
r = [100 92 88];
c = [1 1/80 1];

[rdms,mags] = run_bem_computation(r,c,dippos);
% assertElementsAlmostEqual(rdms, [0.064093 0.064092 0.13532], 'absolute', 1e-3)
% assertElementsAlmostEqual(mags, [1.0498 1.0498 1.0207], 'absolute', 1e-3)
assert(norm(rdms-[0.064093 0.064092 0.13532])<thr)
assert(norm(mags-[1.0498 1.0498 1.0207])<thr)

% 2 Layers, defined from outside to inside
r = [100 92];
c = [1/4 1];

[rdms,mags] = run_bem_computation(r,c,dippos);
% assertElementsAlmostEqual(rdms, [0.15514 0.15514 0.1212], 'absolute', 1e-3)
% assertElementsAlmostEqual(mags, [1.8211 1.8211 1.3606], 'absolute', 1e-3)
assert(norm(rdms-[0.15514 0.15514 0.1212])<thr)
assert(norm(mags-[1.8211 1.8211 1.3606])<thr)

% 1 Layer
r = [100];
c = [1];

[rdms,mags] = run_bem_computation(r,c,dippos);
% assertElementsAlmostEqual(rdms, [0.18934 0.18931 0.0778], 'absolute', 1e-3)
% assertElementsAlmostEqual(mags, [1.3584 1.3583 1.2138], 'absolute', 1e-3)
assert(norm(rdms-[0.18934 0.18931 0.0778])<thr)
assert(norm(mags-[1.3584 1.3583 1.2138])<thr)


  function [rdms,mags] = run_bem_computation(r,c,dippos)
    
    %% Description of the spherical mesh
    [pos, tri] = mesh_sphere(42);
    % [pos, tri] = mesh_sphere(162);
    % [pos, tri] = mesh_sphere(642);
    
    %% Create a set of electrodes on the outer surface
    sens.elecpos = max(r) * pos;
    sens.label = {};
    sens.unit = 'mm';
    nsens = size(sens.elecpos,1);
    for ii=1:nsens
      sens.label{ii} = sprintf('vertex%03d', ii);
    end
    
    %% Create a triangulated mesh, the first boundary is outside
    bnd = [];
    for ii=1:length(r)
      bnd(ii).pos = pos * r(ii);
      bnd(ii).tri = tri;
    end
    
    
    %% Compute the BEM model
    cfg = [];
    cfg.method = 'openmeeg';
    
    cfg.conductivity = c;
    vol_bem = ft_prepare_headmodel(cfg, bnd);
    
    cfg.headmodel = vol_bem;
    cfg.sourcemodel.pos = dippos;
    cfg.sourcemodel.unit = 'mm';
    cfg.elec = sens;
    [grid,cfg] = ft_prepare_leadfield(cfg);
    
    lf_openmeeg = sourcemodel.leadfield{1};
    
    % Rq : ft_compute_leadfield centers the forward fields by default
    % (average reference)
    % lf_openmeeg = lf_openmeeg - repmat(mean(lf_openmeeg),size(lf_openmeeg,1),1);
    
    %% Compute the analytic leadfield
    vol_sphere = [];
    vol_sphere.r = r;
    vol_sphere.cond = c;
    lf_sphere = ft_compute_leadfield(dippos, sens, vol_sphere);
    
    %% Evaluate the quality of the result using RDM and MAG
    rdms = zeros(1,size(lf_openmeeg,2));
    for ii=1:size(lf_openmeeg,2)
      rdms(ii) = norm(lf_openmeeg(:,ii)/norm(lf_openmeeg(:,ii)) - lf_sphere(:,ii) / norm(lf_sphere(:,ii)));
    end
    mags = sqrt(sum(lf_openmeeg.^2))./sqrt(sum(lf_sphere.^2));
    disp(['RDMs: ',num2str(rdms)]);
    disp(['MAGs: ',num2str(mags)]);
    
  end %  function
