function testOpenMEEGeeg
% Test the computation of an EEG leadfield with OpenMEEG

addpath(cd) % Make sure current folder is in the path

%% Set the position of the probe dipole
pos = [0 0 70];

%% Set the radius and conductivities of each of the compartments

% 4 Layers
r = [85 88 92 100];
c = [1 1/20 1/80 1];

[rdms,mags] = run_bem_computation(r,c,pos);
assertElementsAlmostEqual(rdms, [0.019963 0.019962 0.10754], 'absolute', 1e-4)
assertElementsAlmostEqual(mags, [0.84467 0.84469 0.83887], 'absolute', 1e-4)

% 3 Layers
r = [88 92 100];
c = [1 1/80 1];

[rdms,mags] = run_bem_computation(r,c,pos);
assertElementsAlmostEqual(rdms, [0.064093 0.064092 0.13532], 'absolute', 1e-4)
assertElementsAlmostEqual(mags, [1.0498 1.0498 1.0207], 'absolute', 1e-4)

% 2 Layers
r = [92 100];
c = [1 1/4];

[rdms,mags] = run_bem_computation(r,c,pos);
assertElementsAlmostEqual(rdms, [0.15514 0.15514 0.12858], 'absolute', 1e-4)
assertElementsAlmostEqual(mags, [1.8211 1.8211 1.3593], 'absolute', 1e-4)

% 1 Layers
r = [100];
c = [1];

[rdms,mags] = run_bem_computation(r,c,pos);
assertElementsAlmostEqual(rdms, [0.18934 0.18931 0.086394], 'absolute', 1e-4)
assertElementsAlmostEqual(mags, [1.3584 1.3583 1.2138], 'absolute', 1e-4)

end %  function

function [rdms,mags] = run_bem_computation(r,c,pos)

    %% Description of the spherical mesh
    [pnt, tri] = icosahedron42;
    % [pnt, tri] = icosahedron162;
    % [pnt, tri] = icosahedron642;

    %% Create a set of electrodes on the outer surface
    sens.pnt = max(r) * pnt;
    sens.label = {};
    nsens = size(sens.pnt,1);
    for ii=1:nsens
        sens.label{ii} = sprintf('vertex%03d', ii);
    end

    %% Create a BEM volume conduction model
    vol = [];
    for ii=1:length(r)
        vol.bnd(ii).pnt = pnt * r(ii);
        vol.bnd(ii).tri = tri;
    end
    vol.cond = c;

    %% Compute the BEM

    cfg.method = 'openmeeg';

    vol = ft_prepare_bemmodel(cfg, vol);

    cfg.vol = vol;
    cfg.grid.pos = pos;
    cfg.elec = sens;
    grid = ft_prepare_leadfield(cfg);

    lf_openmeeg = grid.leadfield{1};

    % Rq : ft_compute_leadfield centers the forward fields by default
    % (average reference)
    % lf_openmeeg = lf_openmeeg - repmat(mean(lf_openmeeg),size(lf_openmeeg,1),1);

    %% Compute the analytic leadfield
    vol_sphere.r = r;
    vol_sphere.c = c;
    lf_sphere = ft_compute_leadfield(pos, sens, vol_sphere);

    %% Evaluate the quality of the result using RDM and MAG
    rdms = zeros(1,size(lf_openmeeg,2));
    for ii=1:size(lf_openmeeg,2)
        rdms(ii) = norm(lf_openmeeg(:,ii)/norm(lf_openmeeg(:,ii)) - lf_sphere(:,ii) / norm(lf_sphere(:,ii)));
    end
    mags = sqrt(sum(lf_openmeeg.^2))./sqrt(sum(lf_sphere.^2));
    disp(['RDMs: ',num2str(rdms)]);
    disp(['MAGs: ',num2str(mags)]);

end %  function

