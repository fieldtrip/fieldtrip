function fiff_write_dig_file(filename,lpa,nas,rpa,hpi,eeg,eegref,extra)
%
% function fiff_write_dig_file(filename,lpa,nas,rpa,hpi,eeg,eegref,extra)
%
% Create a fif file containing the Polhemus data
%
% filename      Output file name
%
% The following need to be specified in the Neuromag MEG
% coordinate system in meters
%
% lpa           Left auricular point
% nas           Nasion
% rpa           Right auricular point
% hpi           HPI coil locations               (optional)
% eeg           EEG electrode locations          (optional)
% eegref        EEG reference electrode location (optional)
% extra         Additional head surface points   (optional)
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.1  2008/06/04 16:23:24  msh
%   Added fiff_write_dig_file.m
%
%


%
me   = 'MNE:mne_create_dig_file';
%
global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end
%
%   Start writing...
%
fid  = fiff_start_file(filename);
%
% Make the transform to bring LPA and RPA to the x axis and nasion
% to the y axis
%
d1 = nas - lpa;
d2 = rpa - lpa;

alpha = d1*d2'/(d2*d2');
r0 = (1-alpha)*lpa + alpha*rpa;     % Origin is where nasion projects on the LPA - RPA line
ex = d2/sqrt(d2*d2');               % x axis through the auricular points, positive from left to right
ey = (nas-r0);                      % y axis points from the origin to the nasion
ey = ey/sqrt(ey*ey');               % Normalize
ez = cross(ex,ey);                  % Use cross product to get a right-handed coordinate system

T = inv([ ex 0 ; ey 0 ; ez 0 ; r0 1 ]);
% The above gives the coordinate transform from the Neuromag coordinates to the
% Polhemus data coordinates; inverse is needed to go the other way
T = T(:,1:3);   % We just need the first three columns because the result is not an augmented vector

%
% Ready to start
%
fiff_start_block(fid,FIFF.FIFFB_ISOTRAK);

fiff_write_int(fid,FIFF.FIFF_MNE_COORD_FRAME,FIFF.FIFFV_COORD_HEAD);

d.kind  = FIFF.FIFFV_POINT_CARDINAL;
d.ident = FIFF.FIFFV_POINT_LPA;
d.r     = [lpa 1]*T;
fiff_write_dig_point(fid,d);

d.kind  = FIFF.FIFFV_POINT_CARDINAL;
d.ident = FIFF.FIFFV_POINT_NASION;
d.r     = [nas 1]*T;
fiff_write_dig_point(fid,d);

d.kind  = FIFF.FIFFV_POINT_CARDINAL;
d.ident = FIFF.FIFFV_POINT_RPA;
d.r     = [rpa 1]*T;
fiff_write_dig_point(fid,d);


fprintf(1,'Wrote the fiducial locations\n');

if exist('hpi','var')
   for k = 1:size(hpi,1)
      d.kind  = FIFF.FIFFV_POINT_HPI;
      d.ident = k;
      d.r     = [hpi(k,:) 1]*T;
      fiff_write_dig_point(fid,d);
   end
   if (size(hpi,1) > 0)
      fprintf(1,'Wrote %d HPI coil locations\n',size(hpi,1));
   end
end

if exist('eeg','var')
   if exist('eegref','var') && ~isempty(eegref)
      d.kind  = FIFF.FIFFV_POINT_EEG;
      d.ident = 0;
      d.r     = [eegref 1]*T;
      fiff_write_dig_point(fid,d);
      fprintf(1,'Wrote EEG reference electrode location\n');
   end
   for k = 1:size(eeg,1)
      d.kind  = FIFF.FIFFV_POINT_EEG;
      d.ident = k;
      d.r     = [eeg(k,:) 1]*T;
      fiff_write_dig_point(fid,d);
   end
   if (size(eeg,1) > 0)
      fprintf(1,'Wrote %d EEG electrode locations\n',size(eeg,1));
   end
end

if exist('extra','var')
   for k = 1:size(extra,1)
      d.kind  = FIFF.FIFFV_POINT_EXTRA;
      d.ident = k;
      d.r     = [extra(k,:) 1]*T;
      fiff_write_dig_point(fid,d);
   end
   if (size(extra,1) > 0)
      fprintf(1,'Wrote %d extra points\n',size(extra,1));
   end
end

fiff_end_block(fid,FIFF.FIFFB_ISOTRAK);

fiff_end_file(fid);

