function cond = fns_contable_write(varargin)
% Creates the default conductivity table
%
% Use as
%   cond = fns_contable_write
%
% The FNS convention for tissue types is the following:
%     0     "Clear Label"
%     1     "CSF"
%     2     "Gray Matter"
%     3     "White Matter"
%     4     "Fat"
%     5     "Muscle"
%     6     "Muscle/Skin"
%     7     "Skull"
%     8     "Vessels"
%     9     "Around Fat"
%    10     "Dura Matter"
%    11     "Bone Marrow"
%    12     "Eyes"
%

% Copyright (C) 2011, Hung Dang, Cristiano Micheli
tissue     = ft_getopt(varargin, 'tissue', []);
tissueval  = ft_getopt(varargin, 'tissueval', []);
tissuecond = ft_getopt(varargin, 'tissuecond', []);

if isempty(tissue) && isempty(tissueval) && isempty(tissuecond)
  cond = createCondMat;
elseif isempty(tissueval) || isempty(tissuecond)
  error('Both tissue value and conductivity inputs are necessary')
else
  cond = zeros(9,numel(tissueval));
  for i=1:numel(tissueval)
    cond([1,5,9],i) = tissuecond(i);
  end
end

function cond = createCondMat
cond = zeros(9,13);
% Brain tissues index
HH_OUTSIDE        = 0;
HH_CSF            = 1;
HH_GRAY_MATTER    = 2;
HH_WHITE_MATTER   = 3;
HH_FAT            = 4;
HH_MUSCLE         = 5;
HH_SKIN           = 6;
HH_SKULL          = 7;
HH_VESSELS        = 8;
HH_AROUND_FAT     = 9;
HH_DURA           = 10;
HH_BONE_MARROW    = 11;
HH_EYES	          = 12;

% Brain tissue conductivities
HH_OUTSIDE_CON      = 0.00;
HH_CSF_CON          = 1.79;
HH_WHITE_MATTER_CON = 0.14;
HH_GRAY_MATTER_CON  = 0.33;
HH_FAT_CON          = 0.04;
HH_MUSCLE_CON       = 0.11;
HH_SKIN_CON         = 0.44;
HH_SKULL_CON        = 0.018;
HH_VESSELS_CON      = 0.68;
HH_AROUND_FAT_CON   = 0.22;
HH_DURA_CON         = 0.17;
HH_BONE_MARROW_CON  = 0.085;
HH_EYES_CON         = 0.5;

% $$$ Create the conductivity table
cond([1,5,9],HH_OUTSIDE + 1)      = HH_OUTSIDE_CON;
cond([1,5,9],HH_CSF + 1)          = HH_CSF_CON;
cond([1,5,9],HH_WHITE_MATTER + 1) = HH_WHITE_MATTER_CON;
cond([1,5,9],HH_GRAY_MATTER + 1)  = HH_GRAY_MATTER_CON;
cond([1,5,9],HH_FAT + 1)          = HH_FAT_CON;
cond([1,5,9],HH_MUSCLE + 1)       = HH_MUSCLE_CON;
cond([1,5,9],HH_SKULL + 1)        = HH_SKULL_CON;
cond([1,5,9],HH_SKIN + 1)         = HH_SKIN_CON;
cond([1,5,9],HH_VESSELS + 1)      = HH_VESSELS_CON;
cond([1,5,9],HH_AROUND_FAT + 1)   = HH_AROUND_FAT_CON;
cond([1,5,9],HH_DURA + 1)         = HH_DURA_CON;
cond([1,5,9],HH_BONE_MARROW + 1)  = HH_BONE_MARROW_CON;
cond([1,5,9],HH_EYES + 1)         = HH_EYES_CON;
