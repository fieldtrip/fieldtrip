function test_datatype_source

% MEM 1500mb
% WALLTIME 00:10:00

% this function defines a bunch of 'ideal' source structures

% create a set of 3D grid positions
[x,y,z] = ndgrid(-2:5,-3:3,-1:2);
pos     = [x(:) y(:) z(:)];
npos    = size(pos,1);
dim     = size(x);
inside  = false(dim);
inside(2:end-1,2:end-1,:) = true;
outside = find(~inside);
inside  = find(inside);

clear x y z

% create a triangulated mesh
[pnt,tri] = icosahedron162;
npnt      = size(pnt,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the definitions start here 

% 3D grid average time domain
source1     = [];
source1.pos = pos;
source1.dim = dim;
source1.time = (-10:89)./100;
source1.inside  = inside(:)';
source1.outside = outside(:)';
source1.mom     = cell(1,npos);
for k = 1:numel(inside)
  source1.mom{inside(k)} = randn(3,100);
end
source1.momdimord = '{pos}_dir_time';
source1.pow       = source1.mom;
for k = 1:numel(inside)
  source1.pow{inside(k)} = sum(abs(source1.mom{inside(k)}).^2);
end
source1.powdimord = '{pos}_time';

% 3D grid multiple trials time domain
source2 = [];
source2.pos = pos;
source2.dim = dim;
source2.time = (-10:89)./100;
source2.inside  = inside(:)';
source2.outside = outside(:)';
source2.mom     = cell(1,npos);
for k = 1:numel(inside)
  source2.mom{inside(k)} = randn(20,3,100);
end
source2.momdimord = '{pos}_rpt_dir_time';
source2.pow       = source2.mom;
for k = 1:numel(inside)
  source2.pow{inside(k)} = squeeze(sum(abs(source2.mom{inside(k)}).^2,2));
end
source2.powdimord = '{pos}_rpt_time';

% 3D grid multiple trials time domain fixed orientation
source3 = [];
source3.pos = pos;
source3.ori = nan(npos,3); %FIXME or should this be 3xnpos?
source3.ori(inside,:) = randn(numel(inside),3);
for k = 1:numel(inside)
  source3.ori(inside(k),:) = source3.ori(inside(k),:)./norm(source3.ori(inside(k),:));
end
source3.dim = dim;
source3.time = (-10:89)./100;
source3.inside  = inside(:)';
source3.outside = outside(:)';
source3.mom     = cell(1,npos);
for k = 1:numel(inside)
  source3.mom{inside(k)} = randn(20,100);
end
source3.momdimord = '{pos}_rpt_time';
source3.pow       = source2.mom;
for k = 1:numel(inside)
  source3.pow{inside(k)} = sum(abs(source3.mom{inside(k)}).^2,2);
end
source3.powdimord = '{pos}_rpt_time';

% cortical sheet average time domain
source4 = [];
source4.pos = pnt;
source4.tri = tri;
source4.time = (-10:89)./100;
source4.inside  = 1:npnt;
source4.outside = [];
source4.mom     = cell(1,npnt);
inside = source4.inside;
for k = 1:numel(inside)
  source4.mom{inside(k)} = randn(3,100);
end
source4.momdimord = '{pos}_time';
source4.pow       = zeros(npnt,100);
for k = 1:numel(inside)
  source4.pow(k,:) = sum(abs(source4.mom{inside(k)}).^2);
end
source4.powdimord = 'pos_time';

% cortical sheet multiple trials
source5 = [];
source5.pos = pnt;
source5.tri = tri;
source5.time = (-10:89)./100;
source5.inside  = 1:npnt;
source5.outside = [];
source5.mom     = cell(1,npnt);
inside = source5.inside;
for k = 1:numel(inside)
  source5.mom{inside(k)} = randn(20,3,100);
end
source5.momdimord = '{pos}_rpt_dir_time';
source5.pow       = zeros(npnt,20,100);
for k = 1:numel(inside)
  source5.pow(k,:,:) = squeeze(sum(abs(source5.mom{inside(k)}).^2,2));
end
source5.powdimord = 'pos_rpt_time';

% cortical sheet multiple trials with fixed orientation
source6 = [];
source6.pos = pnt;
source6.tri = tri;
source6.ori = nan(npos,3); %FIXME or should this be 3xnpos?
source6.ori(inside,:) = randn(numel(inside),3);
for k = 1:numel(inside)
  source6.ori(inside(k),:) = source6.ori(inside(k),:)./norm(source6.ori(inside(k),:));
end
source6.time = (-10:89)./100;
source6.inside  = 1:npnt;
source6.outside = [];
source6.mom     = cell(1,npnt);
inside = source6.inside;
for k = 1:numel(inside)
  source6.mom{inside(k)} = randn(20,100);
end
source6.momdimord = '{pos}_rpt_time';
source6.pow       = zeros(npnt,20,100);
for k = 1:numel(inside)
  source6.pow(k,:,:) = abs(source6.mom{inside(k)}).^2;
end
source6.powdimord = 'pos_rpt_time';

% 3D grid average freq domain
source7      = [];
source7.pos  = pos;
source7.dim  = dim;
source7.freq = 10;
source7.inside  = inside(:)';
source7.outside = outside(:)';
source7.csd     = cell(1,npos);
for k = 1:numel(inside)
  tmp = randn(3,1)+randn(3,1).*1i;
  source7.csd{inside(k)} = tmp*tmp';
end
source7.csddimord = '{pos}_dir_dir_freq';
source7.pow       = zeros(npos,1);
for k = 1:numel(inside)
  source7.pow(inside(k)) = trace(source7.csd{inside(k)});
end
source7.powdimord = 'pos_freq';

