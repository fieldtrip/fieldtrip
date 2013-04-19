function p = pred_coxphp(gp, x, y, xt, yt, varargin)
%PRED_COXPHP  Integrates the model from zero to point yt (when time is
%  scaled to interval 0-1)
%
%  Description
%    P = PRED_COXPHP(GP,X,Y,XT,YT)
%    If given 1D vector Y, integrates the model from zero to point yt with
%    respect to time. Return P, the probability that event has happened
%    before time yt. YT is vector of size 1xM(or Mx1) indicating times scaled
%    to same interval as the time in timeprocess. Returns matrix P of size NxM
%    where columns correspond to points in YT and rows correspond to rows in X
%    (e.g. people). In case of 2D Y, Integrate model from starting time YT(:,1)
%    to end time YT(:,2). YT is matrix of size Mx2, indicating starting time
%    and end time for every test point.
%

% Copyright (c) 2012-2013 Ville Tolvanen

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

if iscell(gp)
  % prediction for GP_IA cell array
  nGP = numel(gp);
  pp=zeros(size(xt,1),nGP);
  P_TH=zeros(1,nGP);
  for i1=1:nGP
    % make prediction for each gp in cell array
    Gp=gp{i1};
    P_TH(:,i1) = Gp.ia_weight;
    pp(:,i1)=pred_coxphp(Gp, x, y, xt, yt, varargin{:});
  end
  % combine predictions
  p=sum(bsxfun(@times,pp,P_TH),2);
  return
elseif numel(gp.jitterSigma2)>1
  nmc=size(gp.jitterSigma2,1);
  pp=zeros(size(xt,1),nmc);
  for i1=1:nmc
    Gp = take_nth(gp,i1);
    pp(:,i1)=pred_coxphp(Gp, x, y, xt, yt, varargin{:});
  end
  % combine predictions
  p=mean(pp,2);
  return
end

[Ef1, Ef2, Covf] = pred_coxph(gp,x,y,xt, varargin{:});
nsamps = 10000;
ntime=size(gp.lik.stime,2)-1;
if isfield(gp.lik, 'stratificationVariables')
  ind_str=gp.lik.stratificationVariables;
  ux=unique([x(:,ind_str); xt(:,ind_str)],'rows');
  nu=size(ux,1);
  for i1=1:size(ux,1)
    uind{i1}=find(xt(:,ind_str)==ux(i1,:));
  end
  nf1=ntime*nu;
else
  nf1=ntime;
end
sd=gp.lik.stime(2)-gp.lik.stime(1);
Sigm_tmp=Covf;
Sigm_tmp=(Sigm_tmp+Sigm_tmp')./2;
f_star=mvnrnd([Ef1;Ef2], Sigm_tmp, nsamps);

f1=f_star(:,1:nf1);
f2=f_star(:,(nf1+1):end);

la1=exp(f1);
eta2=exp(f2);

mST=zeros(size(eta2,2),1);
if size(y,2) == 1
  if ~isfield(gp.lik,'stratificationVariables')
    % Integrate from zero to yt
    cumsumtmp=cumsum(la1'*sd)';
    for i1=1:size(eta2,2)
      Stime=exp(-bsxfun(@times,cumsumtmp,eta2(:,i1)));
      mStime=mean(Stime);
      if size(yt,1)==size(y,1) && size(yt,2)==size(y,2)
        % Individual integration limits for each input
        mST(i1,1)=1-mStime(binsgeq(gp.lik.xtime,yt(i1)));
      elseif size(yt,1)==1 && size(yt,2)>=1
        % Multiple, but same, integration limits for all of the inputs
        for i=1:size(yt,2)
          mST(i1,i)=1-mStime(binsgeq(gp.lik.xtime,yt(1,i)));
        end
      else
        error('Size of yt is not equal to size of y or 1xT')
      end
    end
    p = mST;
  else
    % Integrate from zero to yt
    for i2=1:nu
      hb=(la1(:,(i2-1)*ntime+1:i2*ntime)'*sd);
      cumsumtmp=cumsum(hb)';
      for i1=1:length(uind{i2})
        ind=uind{i2}(i1);
        Stime=exp(-bsxfun(@times,cumsumtmp,eta2(:,ind)));
        mStime=mean(Stime);
        if size(yt,1)==size(y,1) && size(yt,2)==size(y,2)
          % Individual integration limits for each input
          mST(ind,1)=1-mStime(binsgeq(gp.lik.xtime,yt(ind)));
        elseif size(yt,1)==1 && size(yt,2)>=1
          % Multiple, but same, integration limits for all of the inputs
          for i=1:size(yt,2)
            mST(ind,i)=1-mStime(binsgeq(gp.lik.xtime,yt(1,i)));
          end
        else
          error('Size of yt is not equal to size of y or 1xT')
        end
      end
    end
    p = mST;
  end
  
else
  if size(y,2) ~= size(yt,2)
    error('size(y,2) ~= size(yt,2)');
  end
  if (any(yt(:,2) > gp.lik.stime))
    error('YT has to be scaled to same interval as the timeprocess');
  end
  % Integrate from yt(:,1) to yt(:,2)
  sb=sum(bsxfun(@gt,yt(:,1),gp.lik.stime),2);
  se=sum(bsxfun(@gt,yt(:,2),gp.lik.stime),2);
  
  if ~isfield(gp.lik, 'stratificationVariables')
    for i1 =1:size(eta2,2)
      if sb(i1) ~= se(i1)
        hb=(la1(:,sb(i1)+1:se(i1)-1)'*sd);
        hb=[((gp.lik.stime(sb(i1)+1)-yt(i1,1)).*la1(:,sb(i1)))'; hb; ((yt(i1,2)-gp.lik.stime(se(i1))).*la1(:,se(i1)))'];
      else
        hb = la1(:, se(i1))'*(yt(i1,2) - yt(i1,1));
      end
      cumsumtmp=[zeros(nsamps, sb(i1)) cumsum(hb)'];
      Stime=exp(-bsxfun(@times,cumsumtmp,eta2(:,i1)));
      mStime=mean(Stime);
      mST(i1,1)=1-mStime(end);
    end
    p = mST;
  else
    for i2=1:nu
      for i1=1:size(uind{i2},1)
        ind=uind{i2}(i1);
        nft=(i2-1)*ntime;
        if sb(ind) ~= se(ind)
          hb=(la1(:,(sb(ind)+1+nft):(se(ind)-1+nft))'*sd);
          hb=[((gp.lik.stime(sb(ind)+1)-yt(ind,1)).*la1(:,sb(ind)+nft))'; hb; ((yt(ind,2)-gp.lik.stime(se(ind))).*la1(:,se(ind)+nft))'];
        else
          hb = la1(:, se(ind)+nft)'*(yt(ind,2) - yt(ind,1));
        end
        cumsumtmp=[zeros(nsamps, sb(ind)) cumsum(hb)'];
        Stime=exp(-bsxfun(@times,cumsumtmp,eta2(:,ind)));
        mStime=mean(Stime);
        mST(ind,1)=1-mStime(end);
      end
    end
    p = mST;
  end
end



end
