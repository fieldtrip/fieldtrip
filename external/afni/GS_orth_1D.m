function [m1gs] = GS_orth_1D(fname, rescale, show, oname)
%
%  [A_orth] = GS_orth_1D(fname, rescale, show, oname);
%
%  A function to perform Gram-Schmidt orthogonalization
%  on 1D files.
%
%  fname: Name of AFNI 1D file
%         If nothing is specified, GUI file selector is used.
%         Note that you can pass a matrix instead of fname
%  rescale: If 0 then output matrix is an orthonormal version of
%           the matrix in fname
%           If 1, each column of the orthonormal matrix is scaled
%           by the norm of that column in the original matrix
%           Default is 1
%  show: If 0 then be quiet
%        If 1 then display matrix columns
%        Default is 1
%  oname: If empty then nothing is done.
%         otherwise, write the output matrix to a .1D file of that name
%
%  A_orth: Gram-Schmidt orthogonalized version of the matrix in fname
%
%  This function is a wrapper for Gram_Schmidt_Process.m
%

if (nargin < 4), oname = ''; end
if (nargin < 3), show = 1; end
if (nargin < 2), rescale = 1; end
if (nargin < 1),
   m1 = Read_1D();
else
   if (~isnumeric(fname)),
      m1 = Read_1D(fname);
   else
      m1 = fname;
      fname = '___matrix';
   end
end



nc = size(m1,2);
msize = [ min(nc,8) ,ceil(nc./8)];


if (show),
   figure(1);
   for (k=1:1:nc),
      subplot(msize(1), msize(2), k); plot (m1(:,k));
      if (k==1), title ('Regressors'); end
      ylabel(sprintf('Regressor #%g',k));
      if (~rem(k,msize(1))) xlabel(sprintf('TR')); end
   end
end

%Do the orthonormalization
m1gs = Gram_Schmidt_Process(m1);
ncs = size(m1gs,2);
msize = [ min(ncs,8) ,ceil(ncs./8)];
%rescale back
if (rescale),
   if (ncs == nc),
      for (k=1:1:nc),
         m1gs(:,k) = m1gs(:,k).*norm(m1(:,k));
      end
   else
      fprintf(2,'Rescaling not done, matrix not of full rank\n');
   end
end
if (show),
   figure(2);
   for (k=1:1:ncs),
      subplot(msize(1), msize(2), k); plot (m1gs(:,k));
      if (k==1),
         if (rescale) title ('Orthogonalized Regressors');
         else title ('Orthonormalized Regressors');end
      end
      ylabel(sprintf('Regressor #%g',k));
      if (~rem(k,msize(1))) xlabel(sprintf('TR')); end
   end
end

if (~isempty(oname)),
   wryte3(m1gs, oname);
end
