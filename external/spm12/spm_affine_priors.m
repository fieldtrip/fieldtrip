function [mu,isig] = spm_affine_priors(typ)
% Distribution of the priors used in affine registration
%
% The parameters for this distribution were derived empirically from 227
% scans, that were matched to the ICBM space.
%_______________________________________________________________________
% Copyright (C) 2003-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_affine_priors.m 4155 2011-01-11 15:22:39Z guillaume $

%% Values can be derived by...
%sn = spm_select(Inf,'.*seg_inv_sn.mat$');
%X  = zeros(size(sn,1),12);
%for i=1:size(sn,1),
%    p  = load(deblank(sn(i,:)));
%    M  = p.VF(1).mat*p.Affine/p.VG(1).mat;
%    J  = M(1:3,1:3);
%    V  = sqrtm(J*J');
%    R  = V\J;
%    lV      =  logm(V);
%    lR      = -logm(R);
%    P       = zeros(12,1);
%    P(1:3)  = M(1:3,4);
%    P(4:6)  = lR([2 3 6]);
%    P(7:12) = lV([1 2 3 5 6 9]);
%    X(i,:)  = P';
%end;
%mu   = mean(X(:,7:12));
%XR   = X(:,7:12) - repmat(mu,[size(X,1),1]);
%isig = inv(XR'*XR/(size(X,1)-1))


switch deblank(lower(typ))

case 'mni', % For registering with MNI templates...
    mu   = [0.0667 0.0039 0.0008 0.0333 0.0071 0.1071]';
    isig = 1e4 * [
        0.0902   -0.0345   -0.0106   -0.0025   -0.0005   -0.0163
       -0.0345    0.7901    0.3883    0.0041   -0.0103   -0.0116
       -0.0106    0.3883    2.2599    0.0113    0.0396   -0.0060
       -0.0025    0.0041    0.0113    0.0925    0.0471   -0.0440
       -0.0005   -0.0103    0.0396    0.0471    0.2964   -0.0062
       -0.0163   -0.0116   -0.0060   -0.0440   -0.0062    0.1144];

case 'imni', % For registering with MNI templates...
    mu   = -[0.0667 0.0039 0.0008 0.0333 0.0071 0.1071]';
    isig = 1e4 * [
        0.0902   -0.0345   -0.0106   -0.0025   -0.0005   -0.0163
       -0.0345    0.7901    0.3883    0.0041   -0.0103   -0.0116
       -0.0106    0.3883    2.2599    0.0113    0.0396   -0.0060
       -0.0025    0.0041    0.0113    0.0925    0.0471   -0.0440
       -0.0005   -0.0103    0.0396    0.0471    0.2964   -0.0062
       -0.0163   -0.0116   -0.0060   -0.0440   -0.0062    0.1144];

case 'rigid', % Constrained to be almost rigid...
    mu   = zeros(6,1);
    isig = eye(6)*1e8; % spm_affreg used 1e9

case 'subj', % For inter-subject registration...
    mu   = zeros(6,1);
    isig = 1e3 * [
        0.8876    0.0784    0.0784   -0.1749    0.0784   -0.1749
        0.0784    5.3894    0.2655    0.0784    0.2655    0.0784
        0.0784    0.2655    5.3894    0.0784    0.2655    0.0784
       -0.1749    0.0784    0.0784    0.8876    0.0784   -0.1749
        0.0784    0.2655    0.2655    0.0784    5.3894    0.0784
       -0.1749    0.0784    0.0784   -0.1749    0.0784    0.8876];

case 'eastern', % For East Asian brains to MNI...
    mu   = [0.0719   -0.0040   -0.0032    0.1416    0.0601    0.2578]';
    isig = 1e4 * [
        0.0757    0.0220   -0.0224   -0.0049    0.0304   -0.0327
        0.0220    0.3125   -0.1555    0.0280   -0.0012   -0.0284
       -0.0224   -0.1555    1.9727    0.0196   -0.0019    0.0122
       -0.0049    0.0280    0.0196    0.0576   -0.0282   -0.0200
        0.0304   -0.0012   -0.0019   -0.0282    0.2128   -0.0275
       -0.0327   -0.0284    0.0122   -0.0200   -0.0275    0.0511];

case 'none', % No regularisation...
    mu   = zeros(6,1);
    isig = zeros(6);

otherwise
    error(['"' typ '" not recognised as type of regularisation.']);
end
