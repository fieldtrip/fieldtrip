function [grid] = sel50p(cfg, grid, sens)

% This function will add the field "subspace" to the grid definition.
%
% The subspace projection corresponds to selecting 50% of the 
% channels that are the closest to the dipole.

% set the defaults
if ~isfield(cfg, 'feedback'),   cfg.feedback = 'text';     end

Ninside  = length(sourcemodel.inside);
Nchans   = size(sourcemodel.leadfield{sourcemodel.inside(1)}, 1);

if isfield(sens, 'tra')
  % there is a transformation matrix for magnetometers fo gradiometers
  % weigh each magnetometer coil position with its absolute contribution to
  % the field
  tra = abs(sens.tra);
  sumx = tra * sens.coilpos(:,1);
  sumy = tra * sens.coilpos(:,2);
  sumz = tra * sens.coilpos(:,3);
  pnt = [sumx./sum(tra,2) sumy./sum(tra,2) sumz./sum(tra,2)];
else
  pnt = sens.coilpos;
end

% make a selection matrix for all channels
e = sparse(eye(Nchans));

progress('init', cfg.feedback, 'computing channel selection');
for dipindx=1:Ninside
  % renumber the loop-index variable to make it easier to print the progress bar
  i = sourcemodel.inside(dipindx);
  
  % compute the distance from this dipole to each sensor
  dist = sqrt(sum((pnt-repmat(sourcemodel.pos(i,:), [Nchans 1])).^2, 2));
  
  % define the channels of interest for this dipole
  [dum, indx] = sort(dist);
  sel  = indx(1:round(Nchans/2));
  Nsel = length(sel);
  progress(dipindx/Ninside, 'computing channel selection %d/%d, Nsel=%d\n', dipindx, Ninside, Nsel);
  
  % make a slelection matrix for the 50% nearest-by channels
  grid.subspace{sourcemodel.inside(dipindx)} = e(sel,:);
end
progress('close');

% fill the positions outside the brain with NaNs
for dipindx=sourcemodel.outside(:)'
  grid.subspace{dipindx} = nan;
end
