function test_ft_affinecoordinates

% WALLTIME 01:00:00
% MEM 4gb
% DEPENDENCY ft_affinecoordinates

generic = {
  'als'; 'ali'; 'ars'; 'ari';...
  'pls'; 'pli'; 'prs'; 'pri';...
  'las'; 'lai'; 'ras'; 'rai';...
  'lps'; 'lpi'; 'rps'; 'rpi';...
  'asl'; 'ail'; 'asr'; 'air';...
  'psl'; 'pil'; 'psr'; 'pir';...
  'sal'; 'ial'; 'sar'; 'iar';...
  'spl'; 'ipl'; 'spr'; 'ipr';...
  'sla'; 'ila'; 'sra'; 'ira';...
  'slp'; 'ilp'; 'srp'; 'irp';...
  'lsa'; 'lia'; 'rsa'; 'ria';...
  'lsp'; 'lip'; 'rsp'; 'rip'}';

% the specific ones specify the axes and also the origin
specific = {'ctf', 'bti', 'fourd', 'yokogawa', 'eeglab', 'neuromag', 'itab', 'acpc', 'spm', 'mni', 'fsaverage', 'tal', 'scanras', 'scanlps', 'dicom'};

coordsys = [generic specific];
ncoordsys = numel(coordsys);

f = randi(ncoordsys, 100, 1);
t = randi(ncoordsys, 100, 1);

T1 = zeros(4,4,100);
T2 = zeros(4,4,100);
Tt1 = zeros(100,1);
TT2 = zeros(100,1);
for k = 1:numel(f)
  from = coordsys{f(k)};
  to = coordsys{t(k)};
  fprintf('%s -> %s\n', from, to);
  fprintf('k=%d, new, ', k);
  try
  tic;T1(:,:,k) = ft_affinecoordinates(from,to); Tt1(k,1) = toc;
  fprintf('time elapsed: %2.2g s\n', Tt1(k));
  catch
    fprintf('failed\n');
  end
  fprintf('k=%d, orig, ', k);
  try
  tic;T2(:,:,k) = ft_affinecoordinates_orig(from,to); Tt2(k,1) = toc;
  fprintf('time elapsed: %2.2g s\n', Tt2(k));
  catch
    fprintf('failed\n');
  assert(isequal(T1(:,:,k),T2(:,:,k)));
  end
end


function [transform] = ft_affinecoordinates_orig(original, target, varargin)

% this is the original implementation from Robert

% these are the 48 generic axis orientation triplets, these specify the axes but no origin
%   a = anterior
%   p = posterior
%   l = left
%   r = right
%   s = superior
%   i = inferior

% get the options
requireorigin = ft_getopt(varargin, 'requireorigin', true);

% ensure these are in lower case
original = lower(original);
target = lower(target);

generic = {
  'als'; 'ali'; 'ars'; 'ari';...
  'pls'; 'pli'; 'prs'; 'pri';...
  'las'; 'lai'; 'ras'; 'rai';...
  'lps'; 'lpi'; 'rps'; 'rpi';...
  'asl'; 'ail'; 'asr'; 'air';...
  'psl'; 'pil'; 'psr'; 'pir';...
  'sal'; 'ial'; 'sar'; 'iar';...
  'spl'; 'ipl'; 'spr'; 'ipr';...
  'sla'; 'ila'; 'sra'; 'ira';...
  'slp'; 'ilp'; 'srp'; 'irp';...
  'lsa'; 'lia'; 'rsa'; 'ria';...
  'lsp'; 'lip'; 'rsp'; 'rip'}';

% the specific ones specify the axes and also the origin
specific = {'ctf', 'bti', 'fourd', 'yokogawa', 'eeglab', 'neuromag', 'itab', 'acpc', 'spm', 'mni', 'fsaverage', 'tal', 'scanras', 'scanlps', 'dicom'};

%--------------------------------------------------------------------------
% These approximate alignments are based on transformation matrices that were
% determined by clicking on the CTF fiducial locations in the canonical T1 template
% MRI.
%
% All of the transformation matrices here are expressed in millimeter.

% this is based on the ear canals, see ALIGN_CTF2ACPC
acpc2ctf = [
   0.0000  0.9987  0.0517  34.7467
  -1.0000  0.0000  0.0000   0.0000
   0.0000 -0.0517  0.9987  52.2749
   0.0000  0.0000  0.0000   1.0000
  ];

% this is based on the ear canals, see ALIGN_NEUROMAG2ACPC
acpc2neuromag = [
  1.0000  0.0000  0.0000   0.0000
  0.0000  0.9987  0.0517  34.7467
  0.0000 -0.0517  0.9987  52.2749
  0.0000  0.0000  0.0000   1.0000
  ];

% see http://freesurfer.net/fswiki/CoordinateSystems
fsaverage2mni = [
   0.9975   -0.0073    0.0176   -0.0429
   0.0146    1.0009   -0.0024    1.5496
  -0.0130   -0.0093    0.9971    1.1840
   0.0000    0.0000    0.0000    1.0000
  ];

% this is a 90 degree rotation around the z-axis
ctf2neuromag = [
  0.0000   -1.0000    0.0000    0.0000
  1.0000    0.0000    0.0000    0.0000
  0.0000    0.0000    1.0000    0.0000
  0.0000    0.0000    0.0000    1.0000
  ];

% these are all combinations of 90 degree rotations and/or flips along one of the axes
for i=1:length(generic)
  for j=1:length(generic)
    xxx = generic{i};
    yyy = generic{j};
    eval(sprintf('%s2%s = transform_generic(''%s'', ''%s'');', xxx, yyy, xxx, yyy));
  end
end

% affine transformation from MNI to Talairach, see http://imaging.mrc-cbu.cam.ac.uk/imaging/MniTalairach
% the non-linear (i.e. piecewise linear) transform between MNI and Talairach are implemented elsewhere, see the functions MNI2TAL and TAL2MNI
mni2tal = [
  0.8800    0.0000    0.0000   -0.8000
  0.0000    0.9700    0.0000   -3.3200
  0.0000    0.0500    0.8800   -0.4400
  0.0000    0.0000    0.0000    1.0000
  ];

% the CTF and BTI coordinate system are the same, see http://www.fieldtriptoolbox.org/faq/coordsys/
ctf2bti = eye(4);

% the CTF and EEGLAB coordinate system are the same, see http://www.fieldtriptoolbox.org/faq/coordsys/
ctf2eeglab = eye(4);

% the NEUROMAG and ITAB coordinate system are the same, see http://www.fieldtriptoolbox.org/faq/coordsys/
neuromag2itab = eye(4);

% BTI and 4D are different names for exactly the same system, see http://www.fieldtriptoolbox.org/faq/coordsys/
bti2fourd = eye(4);

% the YOKOGAWA system expresses positions relative to the dewar, not relative to the head
% see https://www.fieldtriptoolbox.org/faq/coordsys/#details-of-the-yokogawa-coordinate-system
yokogawa2als = eye(4);

% the SPM and MNI coordinate system are the same, see http://www.fieldtriptoolbox.org/faq/acpc/
spm2mni = eye(4);

% the SPM (aka MNI) and ACPC coordinate system are not the same, but similar enough, see http://www.fieldtriptoolbox.org/faq/acpc/
spm2acpc = eye(4);
mni2acpc = eye(4);

% the CTF, BTI and 4D coordinate systems are all ALS coordinate systems
% but the origin is poorly defined in ALS, hence converting from ALS to another is problematic
ctf2als   = eye(4);
bti2als   = eye(4);
fourd2als = eye(4);

% the NEUROMAG, ITAB, ACPC, MNI, SPM and FSAVERAGE coordinate systems are all RAS coordinate systems
% but the origin is poorly defined in RAS, hence converting from RAS to another is problematic
neuromag2ras  = eye(4);
itab2ras      = eye(4);
acpc2ras      = eye(4);
mni2ras       = eye(4);
spm2ras       = eye(4);
fsaverage2ras = eye(4);
tal2ras       = eye(4);

% the SCANRAS coordinate system is RAS with the origin at the center opf the gradient coil
scanras2ras     = eye(4);

% the DICOM and SCANLPS coordinate system are the same, and rotated 180 degrees from SCANRAS
dicom2scanlps   = eye(4);
dicom2lps       = eye(4);
scanlps2lps     = eye(4);
scanlps2scanras = lps2ras; % this is a 180 degree rotation around the z-axis

% make the combined and the inverse transformations where possible
coordsys = [specific generic];
implemented = zeros(length(coordsys)); % this is only for debugging
for i=1:numel(coordsys)
  for j=1:numel(coordsys)
    xxx = coordsys{i};
    yyy = coordsys{j};

    if isequal(xxx, yyy)
      % construct the transformations on the diagonal
      eval(sprintf('%s2%s = eye(4);', xxx, yyy));
      implemented(i,j) = 1;
    elseif exist(sprintf('%s2%s', xxx, yyy), 'var')
      % construct the inverse transformations
      eval(sprintf('%s2%s = inv(%s2%s);', yyy, xxx, xxx, yyy));
      implemented(i,j) = 2;
      implemented(j,i) = 2;
    elseif ismember(xxx, specific) && ismember(yyy, generic)
      % try to make the transformation (and inverse) with a two-step approach
      % since we go from specific to generic and thereby loose the origin information anyway, it is fine to use any intermediate step
      for k=1:numel(coordsys)
        zzz = coordsys{k};
        if exist(sprintf('%s2%s', xxx, zzz), 'var') && exist(sprintf('%s2%s', zzz, yyy), 'var')
          eval(sprintf('%s2%s = %s2%s * %s2%s;', xxx, yyy, zzz, yyy, xxx, zzz));
          eval(sprintf('%s2%s = inv(%s2%s);', yyy, xxx, xxx, yyy));
          implemented(i,j) = 3;
          implemented(j,i) = 3;
          break
        end
      end % for k
    elseif ismember(xxx, specific) && ismember(yyy, specific)
      % try to make the transformation (and inverse) with a two-step approach
      % do not use the generic orientation triplets (like RAS and ALS) as intermediate steps between two specific coordinate systems
      for k=1:numel(specific)
        zzz = specific{k};
        if exist(sprintf('%s2%s', xxx, zzz), 'var') && exist(sprintf('%s2%s', zzz, yyy), 'var')
          eval(sprintf('%s2%s = %s2%s * %s2%s;', xxx, yyy, zzz, yyy, xxx, zzz));
          eval(sprintf('%s2%s = inv(%s2%s);', yyy, xxx, xxx, yyy));
          implemented(i,j) = 3;
          implemented(j,i) = 3;
          break
        end
      end % for k
    end

  end % for j
end % for i

% these conversions should be done using FT_VOLUMENORMALISE, as they imply scaling
clear acpc2spm acpc2mni acpc2fsaverage acpc2tal

% converting to/from TAL is only possible for some specific template coordinate systems
clear bti2tal ctf2tal fourd2tal itab2tal neuromag2tal
clear tal2bti tal2ctf tal2fourd tal2itab tal2neuromag

if requireorigin
  % the origin is poorly defined in generic orientation triplets (like RAS and ALS)
  % hence converting between them is problematic with regards to translations
  for i=1:length(generic)
    for j=1:length(specific)
      xxx = generic{i};
      yyy = specific{j};
      eval(sprintf('clear %s2%s', xxx, yyy));
    end
  end
end

if false
  % this is only for checking the coverage of all conversions
  for i=1:length(coordsys)
    for j=1:length(coordsys)
      xxx = coordsys{i};
      yyy = coordsys{j};
      if ~exist(sprintf('%s2%s', xxx, yyy), 'var')
        % update the previous list of implemented transformations, since some have been cleared
        implemented(i,j) = 0;
      end
    end
  end
  figure; imagesc(implemented); caxis([0 3]);
  xticklabels(coordsys); xticks(1:numel(coordsys));
  yticklabels(coordsys); yticks(1:numel(coordsys));
end

if strcmp(original, '4d')
  xxx = 'fourd'; % '4d' is not a valid variable name
else
  xxx = original;
end

if strcmp(target, '4d')
  yyy = 'fourd'; % '4d' is not a valid variable name
else
  yyy = target;
end

if exist(sprintf('%s2%s', xxx, yyy), 'var')
  transform = eval(sprintf('%s2%s', xxx, yyy));
else
  ft_error('converting from %s to %s is not supported', original, target);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to construct generic transformations such as RAS2ALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = transform_generic(from, to)

ap_in  = find(to=='a'   | to=='p');
ap_out = find(from=='a' | from=='p');
lr_in  = find(to=='l'   | to=='r');
lr_out = find(from=='l' | from=='r');
si_in  = find(to=='s'   | to=='i');
si_out = find(from=='s' | from=='i');

% index the axis according to ap,lr,si
order_in  = [ap_in  lr_in  si_in];
order_out = [ap_out lr_out si_out];

% check whether one of the axis needs flipping
flip = 2.*(0.5-double(to(order_in)~=from(order_out)));

T = zeros(4);
for k = 1:3
  T(order_in(k),order_out(k)) = flip(k);
end
T(4,4) = 1;

