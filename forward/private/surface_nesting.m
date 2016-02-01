function order = surface_nesting(bnd, desired)

% SURFACE_NESTING determines what the order of multiple boundaries is to
% get them sorted with the innermost or outermost surface first
%
% Note that it does not check for intersections and may fail for
% intersecting surfaces.

numboundaries = numel(bnd);

% determine the nesting of the compartments
nesting = zeros(numboundaries);
for i=1:numboundaries
  for j=1:numboundaries
    if i~=j
      % determine for a single vertex on each surface if it is inside or outside the other surfaces
      curpos1 = bnd(i).pos(1,:); % any point on the boundary is ok
      curpos  = bnd(j).pos;
      curtri  = bnd(j).tri;
      nesting(i,j) = bounding_mesh(curpos1, curpos, curtri);
    end
  end
end

if sum(nesting(:))~=(numboundaries*(numboundaries-1)/2)
  error('the compartment nesting cannot be determined');
end

if strcmp(desired,'insidefirst')
  % usually the skin will be the outermost, and this should be the first
  % for a three compartment model, the nesting matrix should look like
  %    0 1 1     the first is nested inside the 2nd and 3rd, i.e. the inner skull
  %    0 0 1     the second is nested inside the 3rd, i.e. the outer skull
  %    0 0 0     the third is the most outside, i.e. the skin
  [dum, order] = sort(-sum(nesting,2));
  
elseif strcmp(desired,'outsidefirst')
  % usually the brain (i.e. the inside skull) will be the innermost, and this should be the first
  % for a three compartment model, the nesting matrix should look like
  %    0 0 0     the first is the most outside, i.e. the skin
  %    0 0 1     the second is nested inside the 3rd, i.e. the outer skull
  %    0 1 1     the third is nested inside the 2nd and 3rd, i.e. the inner skull
  [dum, order] = sort(sum(nesting,2));
  
else
  error('unknown surface order "%s"', desired);
end

