function cfg = duneuro_lf(cfg, sens)
% call the post-processing of the LF.
% meg : apply the wight oc sensor (gradiometre, magnetometre)
% eeg : apply the desired reference/ reference == > TODO

% Takfarinas MEDANI, December 2019

if strcmp(cfg.modality,'eeg')
    cfg = postprocess_eeg_lf(cfg); % FIXME THIS NEEDS TO BE FIXED?
end

if strcmp(cfg.modality,'meg')
    cfg = postprocess_meg_lf(cfg, sens);
end

if strcmp(cfg.modality,'meeg')
    cfg = postprocess_meg_lf(cfg, sens);
    cfg = postprocess_eeg_lf(cfg);    
end

function cfg = postprocess_eeg_lf(cfg)

%% substract the mean or not from the electrode
%     cfg.fem_eeg_lf
%     cfg.fem_meg_lf
%TODO : check the minifile parameters and adapt this code
if cfg.lfAvrgRef == 1
    if cfg.useTransferMatrix == 1
        if sum(cfg.lf_fem(1,:)) == 0
            disp('Transforming from elec1 reference to average reference');
            cfg.fem_eeg_lf  = cfg.fem_eeg_lf  - (mean(cfg.fem_eeg_lf,1));
        else
            disp('The average reference is the output of duneuro, please check the mini file');
        end    
    end
end

function cfg = postprocess_meg_lf(cfg, sens)

%% Compute the total magnetic field 
dipoles_pos_ori = [kron(cfg.sourceSpace,ones(3,1)), kron(ones(length(cfg.sourceSpace),1), eye(3))];

% a- Compute the MEG Primary Magnetic field  % apply formula of Sarvas
Bp = compute_B_primary(sens.coilpos, dipoles_pos_ori, sens.coilori);

% b- The total magnetic field B = Bp + Bs;
%  full B-field
Bs =  cfg.meg.Bs;

mu = 4*pi*1e-7; % check the value of the units maybe it needs to be mu = 4*pi*1e-7, I now changed it from 4*pi*1e-4
Bfull = (mu/(4*pi)) * (Bp - Bs);

cfg.meg.lf = Bfull;   

% compute primary magnetic B-field analytically
%
% input:
% coils (Nx3 matrix)
% dipoles (Mx6 matrix)
% projections (Nx3) matrix)

function [Bp] = compute_B_primary(coils, dipoles, projections)

%check input
if size(coils,2)~=3
  error('Column size of coils must be 3.')
end

if size(dipoles,2)~=6
  error('Column size of dipoles must be 6.')
end

if size(projections,2)~=3
  error('Column size of projections must be 3.')
end

% apply formula of Sarvas

dip_pos = dipoles(:,1:3);
dip_mom = dipoles(:,4:6);
Bp = zeros(size(coils,1), size(dipoles,1));
for i = 1:size(coils,1)
  for j = 1 : size(dip_pos,1)
    R = coils(i,:);
    R_0 = dip_pos(j,:);
    A = R - R_0;
    a = norm(A);
    aa = A./(a^3);
    
    BpHelp = cross(dip_mom(j,:),aa);
    Bp(i,j) = BpHelp * projections(i, :)'; % projection of the primary B-field along the coil orientations
  end
end

