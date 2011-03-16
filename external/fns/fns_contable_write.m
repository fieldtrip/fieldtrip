function fns_contable_write(cond,filename)
% Writes the values of conductivity on disk in a csv fashion
% Takes a conductivity matrix as input and if not present defines a standard conductivity matrix
%
% Use as
%   fns_contable_write(cond)
% These values are calculated according to:
% Saleheen HI, Ng KT. New finite difference formulations for general
% inhomogeneous anisotropic bioelectric problems. IEEE Trans Biomed Eng. 1997
% Sep;44(9):800-9

% Copyright (C) 2009, Hung Dang

if isempty(cond) || (size(cond,1)~=9 && size(cond,2)~=12)
  
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
  HH_ELECTRODE      = 20;
  HH_DIPOLE         = 30;
  
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
  cond = zeros(9,12);
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
end

csvwrite(filename,cond);
