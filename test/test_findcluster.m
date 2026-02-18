function test_findcluster

% MEM 6gb
% WALLTIME 00:20:00
% DEPENDENCY findcluster
% DATA no

[dum, ftpath] = ft_version;
cd(fullfile(ftpath, 'private'));

ft_hastoolbox('spm12', 1);

% create a data matrix and a spatial neighbourhoodmatrix
nchan = 5;
C     = diag(ones(1,4),-1)+diag(ones(1,4),1)>0;
dat   = zeros(5,4,3);
dat(1:2,1:2,1:2) = 1;
dat(3,2:3,2:3) = 1;
dat(4,1,1:3) = 1;
dat(5,2:4,3) = 1;
dat = dat>0;

% input is space x something x something
[clus, num] = findcluster(dat,C);
assert(num==3);
assert(nnz(clus(:)==1)==12);
assert(nnz(clus(:)==2)==3);
assert(nnz(clus(:)==3)==3);
clusold = clus;

[clus, num] = findcluster(dat,C,'useoldimplementation',0);
assert(num==3);
assert(nnz(clus(:)==1)==12);
assert(nnz(clus(:)==2)==3);
assert(nnz(clus(:)==3)==3);
assert(isequal(clus,clusold));

% channel 2 is disconnected from the rest
C2 = C;
C2(1,1:3) = [0 0 1];
C2(3,1:3) = [1 0 0];
C2(2,1:3) = 0;
[clus, num] = findcluster(dat,C2);
assert(num==4);
assert(nnz(clus(:)==1)==8);
assert(nnz(clus(:)==2)==4);
assert(nnz(clus(:)==3)==3);
assert(nnz(clus(:)==4)==3);
clusold = clus;

[clus, num] = findcluster(dat,C2,'useoldimplementation',0);
assert(num==4);
assert(nnz(clus(:)==1)==8);
assert(nnz(clus(:)==2)==4);
assert(nnz(clus(:)==3)==3);
assert(nnz(clus(:)==4)==3);
assert(isequal(clus,clusold));

% input is as if it's a spatial 3D reshapeable volume
[clus, num] = findcluster(shiftdim(dat,-1),false);
assert(num==3);
assert(nnz(clus(:)==1)==12);
assert(nnz(clus(:)==2)==3);
assert(nnz(clus(:)==3)==3);
clusold = clus;

[clus, num] = findcluster(shiftdim(dat,-1),false,'useoldimplementation',0);
assert(num==3);
assert(nnz(clus(:)==1)==12);
assert(nnz(clus(:)==2)==3);
assert(nnz(clus(:)==3)==3);
assert(isequal(clus,clusold));

dat4d = repmat(dat, [1 1 1 4]);
dat5d = repmat(dat, [1 1 1 4 5]);
dat6d = repmat(dat, [1 1 1 4 5 6]);
dat7d = repmat(dat, [1 1 1 4 5 6 7]);

[clus4d, num4d] = findcluster(dat4d,C);
assert(num4d==3);
assert(nnz(clus4d(:)==1)==12*4);
assert(nnz(clus4d(:)==2)==3*4);
assert(nnz(clus4d(:)==3)==3*4);
clusold = clus;

[clus4d, num4d] = findcluster(dat4d,C,'useoldimplementation',0);
assert(num4d==3);
assert(nnz(clus4d(:)==1)==12*4);
assert(nnz(clus4d(:)==2)==3*4);
assert(nnz(clus4d(:)==3)==3*4);
assert(isequal(clus,clusold));

[clus5d, num5d] = findcluster(dat5d,C);
assert(num5d==3);
assert(nnz(clus5d(:)==1)==12*4*5);
assert(nnz(clus5d(:)==2)==3*4*5);
assert(nnz(clus5d(:)==3)==3*4*5);
clusold = clus;

[clus5d, num5d] = findcluster(dat5d,C,'useoldimplementation',0);
assert(num5d==3);
assert(nnz(clus5d(:)==1)==12*4*5);
assert(nnz(clus5d(:)==2)==3*4*5);
assert(nnz(clus5d(:)==3)==3*4*5);
assert(isequal(clus,clusold));

[clus6d, num6d] = findcluster(dat6d,C);
assert(num6d==3);
assert(nnz(clus6d(:)==1)==12*4*5*6);
assert(nnz(clus6d(:)==2)==3*4*5*6);
assert(nnz(clus6d(:)==3)==3*4*5*6);
clusold = clus;

[clus6d, num6d] = findcluster(dat6d,C,'useoldimplementation',0);
assert(num6d==3);
assert(nnz(clus6d(:)==1)==12*4*5*6);
assert(nnz(clus6d(:)==2)==3*4*5*6);
assert(nnz(clus6d(:)==3)==3*4*5*6);
assert(isequal(clus,clusold));

[clus7d, num7d] = findcluster(dat7d,C);
assert(num7d==3);
assert(nnz(clus7d(:)==1)==12*4*5*6*7);
assert(nnz(clus7d(:)==2)==3*4*5*6*7);
assert(nnz(clus7d(:)==3)==3*4*5*6*7);
clusold = clus;

[clus7d, num7d] = findcluster(dat7d,C,'useoldimplementation',0);
assert(num7d==3);
assert(nnz(clus7d(:)==1)==12*4*5*6*7);
assert(nnz(clus7d(:)==2)==3*4*5*6*7);
assert(nnz(clus7d(:)==3)==3*4*5*6*7);
assert(isequal(clus,clusold));

[clus4d, num4d] = findcluster(shiftdim(dat4d,-1), false);
assert(num4d==3);
assert(nnz(clus4d(:)==1)==12*4);
assert(nnz(clus4d(:)==2)==3*4);
assert(nnz(clus4d(:)==3)==3*4);
clusold = clus;

[clus4d, num4d] = findcluster(shiftdim(dat4d,-1), false, 'useoldimplementation',0);
assert(num4d==3);
assert(nnz(clus4d(:)==1)==12*4);
assert(nnz(clus4d(:)==2)==3*4);
assert(nnz(clus4d(:)==3)==3*4);
assert(isequal(clus,clusold));

[clus5d, num5d] = findcluster(shiftdim(dat5d,-1), false);
assert(num5d==3);
assert(nnz(clus5d(:)==1)==12*4*5);
assert(nnz(clus5d(:)==2)==3*4*5);
assert(nnz(clus5d(:)==3)==3*4*5);
clusold = clus;

[clus5d, num5d] = findcluster(shiftdim(dat5d,-1), false, 'useoldimplementation',0);
assert(num5d==3);
assert(nnz(clus5d(:)==1)==12*4*5);
assert(nnz(clus5d(:)==2)==3*4*5);
assert(nnz(clus5d(:)==3)==3*4*5);
assert(isequal(clus,clusold));

[clus6d, num6d] = findcluster(shiftdim(dat6d,-1), false);
assert(num6d==3);
assert(nnz(clus6d(:)==1)==12*4*5*6);
assert(nnz(clus6d(:)==2)==3*4*5*6);
assert(nnz(clus6d(:)==3)==3*4*5*6);
clusold = clus;

[clus6d, num6d] = findcluster(shiftdim(dat6d,-1), false, 'useoldimplementation',0);
assert(num6d==3);
assert(nnz(clus6d(:)==1)==12*4*5*6);
assert(nnz(clus6d(:)==2)==3*4*5*6);
assert(nnz(clus6d(:)==3)==3*4*5*6);
assert(isequal(clus,clusold));

[clus7d, num7d] = findcluster(shiftdim(dat7d,-1), false);
assert(num7d==3);
assert(nnz(clus7d(:)==1)==12*4*5*6*7);
assert(nnz(clus7d(:)==2)==3*4*5*6*7);
assert(nnz(clus7d(:)==3)==3*4*5*6*7);
clusold = clus;

[clus7d, num7d] = findcluster(shiftdim(dat7d,-1), false, 'useoldimplementation',0);
assert(num7d==3);
assert(nnz(clus7d(:)==1)==12*4*5*6*7);
assert(nnz(clus7d(:)==2)==3*4*5*6*7);
assert(nnz(clus7d(:)==3)==3*4*5*6*7);
assert(isequal(clus,clusold));


% some comparisons (+timing) for the old and new style implementation
load ctf275_neighb
cfg            = [];
cfg.channel    = {neighbours.label}';
cfg.neighbours = neighbours;
C = channelconnectivity(cfg);

krn = repmat(shiftdim(gausswin(10),-2),[10 10 1]) .* ...
  repmat(shiftdim(gausswin(10),-1),[10 1 10]) .* ...
  repmat(gausswin(10),[1 10 10]);
krn = krn./max(krn(:));


t1 = zeros(10,1);
t2 = zeros(10,1);
deltan = zeros(10,1);
for k = 1:10
  dat = randn(275,50,60);
  datc = convn(dat,krn,'same');
  onoff = datc>15;
  tic;[c1, n1] = findcluster(onoff, C); t1(k) = toc;
  tic;[c2, n2] = findcluster(onoff, C, 'useoldimplementation', 0); t2(k) = toc;
  %assert(isequal(double(c1),c2)); % this is too strict, since the
  %numbering can be different, and - also - the old implementation tends to
  %'miss' to connect isolated pixels in a slice that connects spatially
  %to other 'slices'
  assert(n2<=n1);
  deltan(k) = n1-n2;

end

t1 = zeros(10,1);
t2 = zeros(10,1);
deltan = zeros(10,1);
for k = 1:10
  dat = randn(275,10,12);
  datc = convn(dat,krn,'same');
  onoff = datc>15;
  tic;[c1, n1] = findcluster(onoff, C); t1(k) = toc;
  tic;[c2, n2] = findcluster(onoff, C, 'useoldimplementation', 0); t2(k) = toc;
  %assert(isequal(double(c1),c2)); % this is too strict, since the
  %numbering can be different, and - also - the old implementation tends to
  %'miss' to connect isolated pixels in a slice that connects spatially
  %to other 'slices'
  assert(n2<=n1);
  deltan(k) = n1-n2;
end

ctx = ft_read_headshape(fullfile(ftpath, 'template/sourcemodel', 'cortex_8196.surf.gii'));
C   = triangle2connectivity(ctx.tri);
t1 = zeros(10,1);
t2 = zeros(10,1);
deltan = zeros(10,1);
for k = 1:10
  dat = randn(8196,20);
  datc = convn(dat,krn(:,:,2),'same');
  onoff = datc>1.5;
  tic;[c1, n1] = findcluster(onoff, C); t1(k) = toc;
  tic;[c2, n2] = findcluster(onoff, C, 'useoldimplementation', 0); t2(k) = toc;
  %assert(isequal(double(c1),c2)); % this is too strict, since the
  %numbering can be different, and - also - the old implementation tends to
  %'miss' to connect isolated pixels in a slice that connects spatially
  %to other 'slices'
  assert(n2<=n1);
  deltan(k) = n1-n2;
end

t1 = zeros(10,1);
t2 = zeros(10,1);
ctx = ft_read_headshape(fullfile(ftpath, 'template/sourcemodel', 'cortex_20484.surf.gii'));
C   = triangle2connectivity(ctx.tri);

deltan = zeros(10,1);
for k = 1:10
  dat = randn(20484,20);
  datc = convn(dat,krn(:,:,2),'same');
  onoff = datc+randn(20484,1)*linspace(-1,1,20)>2.5;
  tic;[c1, n1] = findcluster(onoff, C); t1(k) = toc;
  tic;[c2, n2] = findcluster(onoff, C, 'useoldimplementation', 0); t2(k) = toc;
  %assert(isequal(double(c1),c2)); % this is too strict, since the
  %numbering can be different, and - also - the old implementation tends to
  %'miss' to connect isolated pixels in a slice that connects spatially
  %to other 'slices'
  assert(n2<=n1);
  deltan(k) = n1-n2;
end

t1 = zeros(10,1);
t2 = zeros(10,1);

deltan = zeros(10,1);
for k = 1:10
  dat = randn(20484,20,20);
  datc = convn(dat,krn,'same');
  onoff = datc>20;
  tic;[c1, n1] = findcluster(onoff, C); t1(k) = toc;
  tic;[c2, n2] = findcluster(onoff, C, 'useoldimplementation', 0); t2(k) = toc;
  %assert(isequal(double(c1),c2)); % this is too strict, since the
  %numbering can be different, and - also - the old implementation tends to
  %'miss' to connect isolated pixels in a slice that connects spatially
  %to other 'slices'
  assert(n2<=n1);
  deltan(k) = n1-n2;
end


t1 = zeros(10,1);
t2 = zeros(10,1);

deltan = zeros(10,1);
for k = 1:10
  dat = randn(20484,1,1);
  datc = convn(dat,krn(:,5,5), 'same');
  onoff = datc>3;
  tic;[c1, n1] = findcluster(onoff, C); t1(k) = toc;
  tic;[c2, n2] = findcluster(onoff, C, 'useoldimplementation', 0); t2(k) = toc;
  %assert(isequal(double(c1),c2)); % this is too strict, since the
  %numbering can be different, and - also - the old implementation tends to
  %'miss' to connect isolated pixels in a slice that connects spatially
  %to other 'slices'
  assert(n2<=n1);
  deltan(k) = n1-n2;
end







