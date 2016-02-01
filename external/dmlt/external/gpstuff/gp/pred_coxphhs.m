function [mh,mS] = pred_coxphhs(gp, x, y, xt, varargin)
%PRED_COXPHHS  Return hazard and survival functions
%
%  Description
%    [H,S] = PRED_COXPHHS(GP,X,Y,XT) 
%    Returns hazard and survival functions for the inputs XT. GP is the
%    Gaussian process structure, X and Y are the training inputs and
%    outputs.
%

% Copyright (c) 2012-2013 Ville Tolvanen, Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

ntime=size(gp.lik.xtime,1);
if iscell(gp)
  % prediction for GP_IA cell array
  nGP = numel(gp);
  mh=zeros(size(xt,1),ntime,nGP);
  mS=zeros(size(xt,1),ntime,nGP);
  P_TH=zeros(1,nGP);
  for i1=1:nGP
    % make prediction for each gp in cell array
    Gp=gp{i1};
    P_TH(:,i1) = Gp.ia_weight;
    [mh(:,:,i1), mS(:,:,i1)]=pred_coxphhs(Gp, x, y, xt, varargin{:});
    mh(:,:,i1)=P_TH(:,i1)*mh(:,:,i1);
    mS(:,:,i1)=P_TH(:,i1)*mS(:,:,i1);
  end
  % combine predictions
  mh=sum(mh,3);
  mS=sum(mS,3);
  return
elseif numel(gp.jitterSigma2)>1
  nmc=size(gp.jitterSigma2,1);
  mh=zeros(size(xt,1),ntime,nmc);
  mS=zeros(size(xt,1),ntime,nmc);
  for i1=1:nmc
    Gp = take_nth(gp,i1);
    [mh(:,:,i1), mS(:,:,i1)]=pred_coxphhs(Gp, x, y, xt, varargin{:});
  end
  % combine predictions
  mh=mean(mh,3);
  mS=mean(mS,3);
  return
end

[Ef1, Ef2, Covf] = pred_coxph(gp,x,y,xt, varargin{:});
nsamps = 10000;  
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
% f_star=mvnrnd(Ef1, Sigm_tmp(1:ntime,1:ntime), nsamps);
f_star=mvnrnd([Ef1;Ef2], Sigm_tmp, nsamps);

f1=f_star(:,1:nf1);
f2=f_star(:,(nf1+1):end);

la1=exp(f1);
eta2=exp(f2);

if ~isfield(gp.lik, 'stratificationVariables')
  hb=(la1'*sd);
  cumsumtmp=cumsum(hb)';
  for i1=1:size(eta2,2)
    mh(i1,:)=mean(bsxfun(@times,hb',eta2(:,i1)));
    Stime=exp(-bsxfun(@times,cumsumtmp,eta2(:,i1)));
    mS(i1,:)=mean(Stime);
  end
else
  for i2=1:length(uind)
    hb=(la1(:,(i2-1)*ntime+1:i2*ntime)'*sd);
    cumsumtmp=cumsum(hb)';
    for i1=1:size(uind{i2},1)
      ind=uind{i2}(i1);
      mh(ind,:)=mean(bsxfun(@times,hb',eta2(:,ind)));
      Stime=exp(-bsxfun(@times,cumsumtmp,eta2(:,ind)));
      mS(ind,:)=mean(Stime);
    end
  end
end



end
