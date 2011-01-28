function lf=openmeeg_helper(pos,vol,splitdim)
%
% Splits the calculation of leadfields in smaller chunks
% pos: positions of dipoles
% vol: volume conductor
% splitdim: splitting factor 

lf   = [];
ndip = size(pos,1);
numsol   = fix(ndip./splitdim);
nleft    = rem(ndip,splitdim);
lf = zeros(size(vol.mat,1),ndip*3);
if isfield(vol,'mat')
    for i = 1:splitdim
      dsm = openmeeg_dsm(pos((numsol*(i - 1) + 1):numsol*i,:),vol);
      tmp = vol.mat*dsm;
      lf(:,(numsol*(i - 1)*3 + 1):numsol*i*3) = tmp;
    end
    dsm = openmeeg_dsm(pos(ndip-nleft+1:ndip,:),vol);
    tmp = vol.mat*dsm;
    lf(:,(ndip-nleft)*3+1:ndip*3) = tmp;
else
  error('No system matrix is present, BEM head model not calculated yet')
end
