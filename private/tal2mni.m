function outpoints = tal2mni(inpoints)
% Converts coordinates to MNI brain best guess
% from Talairach coordinates
% FORMAT outpoints = tal2mni(inpoints)
% Where inpoints is N by 3 or 3 by N matrix of coordinates
%  (N being the number of points)
% outpoints is the coordinate matrix with MNI points
% Matthew Brett 2/2/01

% ensure that SPM is available, needed for spm_matrix
hasspm = ft_hastoolbox('spm8up', 3) || ft_hastoolbox('spm2', 1);

dimdim = find(size(inpoints) == 3);
if isempty(dimdim)
  ft_error('input must be a N by 3 or 3 by N matrix')
end
if dimdim == 2
  inpoints = inpoints';
end

% Transformation matrices, different zooms above/below AC
rotn  = spm_matrix([0 0 0 0.05]);
upz   = spm_matrix([0 0 0 0 0 0 0.99 0.97 0.92]);
downz = spm_matrix([0 0 0 0 0 0 0.99 0.97 0.84]);

inpoints = [inpoints; ones(1, size(inpoints, 2))];
% Apply inverse translation
inpoints = inv(rotn)*inpoints;

tmp = inpoints(3,:)<0;  % 1 if below AC
inpoints(:, tmp) = inv(downz) * inpoints(:, tmp);
inpoints(:, ~tmp) = inv(upz) * inpoints(:, ~tmp);
outpoints = inpoints(1:3, :);
if dimdim == 2
  outpoints = outpoints';
end




