function [vol] = yokogawa2vol(hdr);

% YOKOGAWA2VOL converts a spherical volume conductor model that can
% be present in the header of a datafile into a structure that can
% be used by FieldTrip.

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: yokogawa2vol.m,v $
% Revision 1.1  2005/09/06 08:54:30  roboos
% new implementations for the Yokogawa 160 channel MEG system
%

% hdr = read_yokogawa_header(filename);
hdr = hdr.orig; % use the original Yokogawa header, not the FieldTrip header

if isfield(hdr.mri_info, 'model_name') && strcmp(hdr.mri_info.model_name, 'SphericalModel')
  % single sphere volume conduction model
  vol.r = hdr.mri_info.r;
  vol.o = [hdr.mri_info.cx hdr.mri_info.cy hdr.mri_info.cz];
else
  error('unsupported volume conductor model');
end
