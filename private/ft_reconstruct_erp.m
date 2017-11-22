function [orig, reconstructed] = ft_reconstruct_erp(input, avgflag)

if nargin<2
  avgflag = true;
end

[nsmp,ncomp]  = size(input.erp_est);
ntrl          = size(input.amp_est,1);
rejectflag    = input.rejectflag;

if avgflag
reconstructed = zeros(nsmp,1);
orig          = zeros(nsmp,1);
for k = 1:ntrl
  if ~rejectflag(k)
    orig = orig + input.data_init(:,k);
    
    for m=1:ncomp
      tmp           = input.erp_est(:,m);
      reconstructed = reconstructed + input.amp_est(k,m).*fun_shift(tmp, input.lat_est(k,m), 1);
    end
  end
end
reconstructed = reconstructed./sum(rejectflag==0);
orig          = orig./sum(rejectflag==0);
else
  reconstructed = zeros(nsmp,sum(rejectflag==0));
  orig          = input.data_init(:,rejectflag==0);
  indx = 0;
  for k = 1:ntrl
    if ~rejectflag(k)
      indx = indx+1;
      for m=1:ncomp
        tmp                   = input.erp_est(:,m);
        reconstructed(:,indx) = reconstructed(:,indx) + input.amp_est(k,m).*fun_shift(tmp, input.lat_est(k,m), 1);
      end
    end
  end
end