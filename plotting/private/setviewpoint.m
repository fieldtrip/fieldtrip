function setviewpoint(ax, coordsys, viewpoint)

% SETVIEWPOINT changes the viewpoint for a 3D image that contains data in a known coordinate system
%
% Use as
%   setviewpoint(ax, coordsys, viewpoint)
%
% For example
%   setviewpoint(gca, 'mni', 'left')
%
% See alo GETORTHOVIEWPOS, COORDSYS2LABEL

if isempty(coordsys)
  coordsys = 'unknown';
end

switch viewpoint
  case 'front'
    viewpoint = 'anterior';
  case 'back'
    viewpoint = 'posterior';
  case 'top'
    viewpoint = 'superior';
  case 'bottom'
    viewpoint = 'inferior';
end

switch lower(coordsys)
  case {'ras' 'scanras' 'nifti' 'neuromag' 'itab' 'acpc' 'spm' 'mni' 'tal'}
    switch viewpoint
      case 'superior'
        view(ax, [0 0 1]); % the nose is pointing up
      case 'inferior'
        view(ax, [180 -90]); % not exactly the same as [0 0 -1], this causes the nose pointing up
      case 'left'
        view(ax, [-1 0 0]);
      case 'right'
        view(ax, [1 0 0]);
      case 'anterior'
        view(ax, [0 1 0]);
      case 'posterior'
        view(ax, [0 -1 0]);
    end % switch viewpoint
  case {'als' 'ctf' '4d' 'bti' 'eeglab'}
    switch viewpoint
      case 'superior'
        view(ax, [-90 90]); % not exactly the same as [0 0 1], this causes the nose pointing up
      case 'inferior'
        view(ax, [90 -90]); % not exactly the same as [0 0 -1], this causes the noise pointing up
      case 'left'
        view(ax, [0 1 0]);
      case 'right'
        view(ax, [0 -1 0]);
      case 'anterior'
        view(ax, [1 0 0]);
      case 'posterior'
        view(ax, [-1 0 0]);
    end % switch viewpoint
  case {'lps' 'scanlps' 'dicom'}
    switch viewpoint
      case 'superior'
        view(ax, [180 90]); % not exactly the same as [0 0 1], this causes the nose pointing up
      case 'inferior'
        view(ax, [0 -90]); % not exactly the same as [0 0 -1], this causes the noise pointing up
      case 'left'
        view(ax, [1 0 0]);
      case 'right'
        view(ax, [-1 0 0]);
      case 'anterior'
        view(ax, [0 -1 0]);
      case 'posterior'
        view(ax, [0 1 0]);
    end % switch viewpoint
  case {'rsp' 'paxinos'}
    switch viewpoint
      case 'superior'
        view(ax, [0 1 0]);
      case 'inferior'
        view(ax, [0 -1 0]);
      case 'left'
        view(ax, [-1 0 0]);
      case 'right'
        view(ax, [1 0 0]);
      case 'anterior'
        view(ax, [0 0 -1]);
      case 'posterior'
        view(ax, [0 0 1]);
    end % switch viewpoint
  otherwise
    ft_warning('cannot change the viewpoint for the "%s" coordinate system', coordsys);
end % switch coordsys
