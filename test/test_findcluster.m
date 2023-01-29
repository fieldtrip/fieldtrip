function test_findcluster

% MEM 8gb
% WALLTIME 00:20:00
% DEPENDENCY findcluster

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
assert(sum(clus(:)==1)==12);
assert(sum(clus(:)==2)==3);
assert(sum(clus(:)==3)==3);

% channel 2 is disconnected from the rest
C2 = C;
C2(1,1:3) = [0 0 1];
C2(3,1:3) = [1 0 0];
C2(2,1:3) = 0;
[clus, num] = findcluster(dat,C2);
assert(num==4);
assert(sum(clus(:)==1)==8);
assert(sum(clus(:)==2)==4);
assert(sum(clus(:)==3)==3);
assert(sum(clus(:)==4)==3);

% input is as if it's a spatial 3D reshapeable volume
[clus, num] = findcluster(shiftdim(dat,-1),false);
assert(num==3);
assert(sum(clus(:)==1)==12);
assert(sum(clus(:)==2)==3);
assert(sum(clus(:)==3)==3);

dat4d = repmat(dat, [1 1 1 4]);
dat5d = repmat(dat, [1 1 1 4 5]);
dat6d = repmat(dat, [1 1 1 4 5 6]);
dat7d = repmat(dat, [1 1 1 4 5 6 7]);

[clus4d, num4d] = findcluster(dat4d,C);
assert(num4d==3);
assert(sum(clus4d(:)==1)==12*4);
assert(sum(clus4d(:)==2)==3*4);
assert(sum(clus4d(:)==3)==3*4);

[clus5d, num5d] = findcluster(dat5d,C);
assert(num5d==3);
assert(sum(clus5d(:)==1)==12*4*5);
assert(sum(clus5d(:)==2)==3*4*5);
assert(sum(clus5d(:)==3)==3*4*5);

[clus6d, num6d] = findcluster(dat6d,C);
assert(num6d==3);
assert(sum(clus6d(:)==1)==12*4*5*6);
assert(sum(clus6d(:)==2)==3*4*5*6);
assert(sum(clus6d(:)==3)==3*4*5*6);

[clus7d, num7d] = findcluster(dat7d,C);
assert(num7d==3);
assert(sum(clus7d(:)==1)==12*4*5*6*7);
assert(sum(clus7d(:)==2)==3*4*5*6*7);
assert(sum(clus7d(:)==3)==3*4*5*6*7);

[clus4d, num4d] = findcluster(shiftdim(dat4d,-1), false);
assert(num4d==3);
assert(sum(clus4d(:)==1)==12*4);
assert(sum(clus4d(:)==2)==3*4);
assert(sum(clus4d(:)==3)==3*4);

[clus5d, num5d] = findcluster(shiftdim(dat5d,-1), false);
assert(num5d==3);
assert(sum(clus5d(:)==1)==12*4*5);
assert(sum(clus5d(:)==2)==3*4*5);
assert(sum(clus5d(:)==3)==3*4*5);

[clus6d, num6d] = findcluster(shiftdim(dat6d,-1), false);
assert(num6d==3);
assert(sum(clus6d(:)==1)==12*4*5*6);
assert(sum(clus6d(:)==2)==3*4*5*6);
assert(sum(clus6d(:)==3)==3*4*5*6);

[clus7d, num7d] = findcluster(shiftdim(dat7d,-1), false);
assert(num7d==3);
assert(sum(clus7d(:)==1)==12*4*5*6*7);
assert(sum(clus7d(:)==2)==3*4*5*6*7);
assert(sum(clus7d(:)==3)==3*4*5*6*7);

