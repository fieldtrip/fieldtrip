function pos = getorthoviewpos(pos, coordsys, viewpoint)

% GETORTHOVIEWPOS obtains the orthographic projections of 3D positions
% based on a given coordinate system and viewpoint
%
% Use as
%   getorthoviewpos(pos, coordsys, viewpoint)
%
% For example
%   getorthoviewpoint(pos, 'mni', 'superior')
%
% See alo SETVIEWPOINT, COORDSYS2LABEL

if size(pos,2)~=3
  ft_error('XYZ coordinates are required to obtain the orthographic projections based on a viewpoint')
end

% create view(az,el) transformation matrix
switch coordsys
  case {'ras' 'neuromag' 'itab' 'acpc' 'spm' 'mni' 'tal'}
    switch viewpoint
      case 'left'
        transmat = viewmtx(-90, 0);
      case 'right'
        transmat = viewmtx(90, 0);
      case 'topleft'
        transmat = viewmtx(-90, 45);
      case 'topright'
        transmat = viewmtx(90, 45);
      case 'superior'
        transmat = viewmtx(0, 90);
      case 'inferior'
        transmat = viewmtx(180, -90);
      case 'posterior'
        transmat = viewmtx(0, 0);
      case 'anterior'
        transmat = viewmtx(180, 0);
      otherwise
        ft_error('orthographic projection using viewpoint "%s" is not supported', viewpoint)
    end % switch viewpoint
  case {'als' 'ctf' '4d' 'bti'}
    switch viewpoint
      case 'left'
        transmat = viewmtx(180, 0);
      case 'right'
        transmat = viewmtx(0, 0);
      case 'topleft'
        transmat = viewmtx(180, 45);
      case 'topright'
        transmat = viewmtx(0, 45);
      case 'superior'
        transmat = viewmtx(-90, 90);
      case 'inferior'
        transmat = viewmtx(90, -90);
      case 'posterior'
        transmat = viewmtx(-90, 0);
      case 'anterior'
        transmat = viewmtx(90, 0);
      otherwise
        ft_error('orthographic projection using viewpoint "%s" is not supported', viewpoint)
    end % switch viewpoint
  otherwise
    ft_error('orthographic projection using coordinate system "%s" is not supported', coordsys)
end % switch coordsys

% extract xy
pos      = ft_warp_apply(transmat, pos, 'homogenous');
pos      = pos(:,[1 2]);