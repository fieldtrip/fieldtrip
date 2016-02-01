function lik = lik_coxph(varargin)
%LIK_COXPH    Create a Cox proportional hazard likelihood structure
%
%  Description
%    LIK = LIK_COXPH('PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    creates a proportional hazard model where a piecewise log-constant
%    baseline hazard is assumed.  
%    
%    The likelihood contribution for the ith observation is
%
%      l_i = h_i(y_i)^(1-z_i)*[exp(-int_0^y_i*h_i dt)],
%    
%    where hazard is h_i=h_0(y_i)*exp(f_i). A zero mean Gaussian process
%    prior is placed for f = [f_1, f_2,...,f_n] ~ N(0, C). C is the
%    covariance matrix, whose elements are given as C_ij = c(x_i, x_j |
%    th). The function c(x_i, x_j| th) is covariance function and th its
%    parameters, hyperparameters. We place a hyperprior for
%    hyperparameters, p(th). 
%
%    The time axis is partioned into K intervals with equal lengths:
%    0 = s_0 < s_1 < ... < s_K, where s_K > y_i for all i. The baseline
%    hazard rate function h_0 is piecewise constant,  
%
%      h_0(t) = la_k,
%
%    when t belongs to the interval (s_{k-1},s_k] and where ft_k=log(la_k).
%    The hazard rate function is smoothed by assuming another Gaussian
%    process prior ft = [ft_1, ft_2,...,ft_K] ~ N(0, C). 
%
%    z is a vector of censoring indicators with z = 0 for uncensored event
%    and z = 1 for right censored event. 
%
%    When using the Coxph likelihood you need to give the vector z
%    as an extra parameter to each function that requires also y. 
%    For example, you should call gpla_e as follows: gpla_e(w, gp,
%    x, y, 'z', z)
%
%  See also
%    GP_SET, LIK_*, PRIOR_*
%

% Copyright (c) 2007-2010 Jarno Vanhatalo & Jouni Hartikainen
% Copyright (c) 2010 Aki Vehtari
% Copyright (c) 2011 Jaakko Riihimäki

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'LIK_COXPH';
  ip.addOptional('lik', [], @isstruct);
  ip.addParamValue('S', linspace(0,1.001,50), @(x) isvector(x));
  ip.addParamValue('stratificationVariables', [], @(x) isvector(x) && all(rem(x,1)==0));
  ip.addParamValue('removeStratificationVariables', [], @(x) ismember(x, {'on' 'off'}));
  ip.parse(varargin{:});
  lik=ip.Results.lik;
  
  if isempty(lik)
    init=true;
    lik.type = 'Coxph';
    lik.nondiagW=true;
  else
    if ~isfield(lik,'type') || ~isequal(lik.type,'Coxph')
      error('First argument does not seem to be a valid likelihood function structure')
    end
    init=false;
  end
  
  if init || ~ismember('S', ip.UsingDefaults)
    s=ip.Results.S;
    xtmp=zeros(length(s)-1,1);
    for i2=1:(length(s)-1)
      xtmp(i2,1)=mean([s(i2) s(i2+1)]);
    end
    lik.xtime=xtmp;
    lik.stime=s;
  end
  
  if ~ismember('stratificationVariables', ip.UsingDefaults)
    lik.stratificationVariables = ip.Results.stratificationVariables;
    if ~ismember('stratificationVariables', ip.UsingDefaults)
      if isequal(ip.Results.removeStratificationVariables, 'on')
        lik.removeStratificationVariables=true;
      end
    end
  end
  
  if init
    % Set the function handles to the nested functions
    lik.fh.pak = @lik_coxph_pak;
    lik.fh.unpak = @lik_coxph_unpak;
    lik.fh.ll = @lik_coxph_ll;
    lik.fh.llg = @lik_coxph_llg;    
    lik.fh.llg2 = @lik_coxph_llg2;
    lik.fh.llg3 = @lik_coxph_llg3;
    lik.fh.tiltedMoments = @lik_coxph_tiltedMoments;
    lik.fh.predy = @lik_coxph_predy;
    lik.fh.invlink = @lik_coxph_invlink;
    lik.fh.recappend = @lik_coxph_recappend;
    lik.fh.predcdf= @lik_coxph_predcdf;
  end

  function [w,s] = lik_coxph_pak(lik)
  %LIK_COXPH_PAK  Combine likelihood parameters into one vector.
  %
  %  Description 
  %    W = LIK_COXPH_PAK(LIK) takes a likelihood structure LIK and
  %    combines the parameters into a single row vector W. This is 
  %    a mandatory subfunction used for example in energy and 
  %    gradient computations.
  %     
  %       w = log(lik.disper)
  %
  %   See also
  %   LIK_COXPH_UNPAK, GP_PAK
    
    w=[];s={};
  end


  function [lik, w] = lik_coxph_unpak(lik, w)
  %LIK_COXPH_UNPAK  Extract likelihood parameters from the vector.
  %
  %  Description
  %    [LIK, W] = LIK_COXPH_UNPAK(W, LIK) takes a likelihood
  %    structure LIK and extracts the parameters from the vector W
  %    to the LIK structure. This is a mandatory subfunction used 
  %    for example in energy and gradient computations.
  %     
  %   Assignment is inverse of  
  %       w = log(lik.disper)
  %
  %   See also
  %   LIK_COXPH_PAK, GP_UNPAK
  
    lik=lik;
    w=w;
    
  end
  
  function ll = lik_coxph_ll(lik, y, f, z)
  %LIK_COXPH_LL  Log likelihood
  %
  %  Description
  %    LL = LIK_COXPH_LL(LIK, Y, F, Z) takes a likelihood
  %    structure LIK, incedence counts Y, expected counts Z, and
  %    latent values F. Returns the log likelihood, log p(y|f,z).
  %    This subfunction is needed when using Laplace approximation
  %    or MCMC for inference with non-Gaussian likelihoods. This 
  %    subfunction is also used in information criteria (DIC, WAIC)
  %    computations.
  %
  %  See also
  %    LIK_COXPH_LLG, LIK_COXPH_LLG3, LIK_COXPH_LLG2, GPLA_E
    
    if isempty(z)
      error(['lik_coxph -> lik_coxph_ll: missing z!    '... 
             'Coxph likelihood needs the expected number of    '...
             'occurrences as an extra input z. See, for         '...
             'example, lik_coxph and gpla_e.               ']);
    end
    
    [n,ny]=size(y);    
    ntime=size(lik.xtime,1);
    i3v=ones(size(y,1),1);
    if isfield(lik, 'stratificationVariables')
      f1=f(1:ntime*lik.n_u);
      f2=f((ntime*lik.n_u+1):(ntime*lik.n_u+n));
      for i=1:lik.n_u
        i3v(lik.stratind{i})=i;
      end
    else
      f1=f(1:ntime);
      f2=f((ntime+1):(ntime+n));
    end    
    
    la1=exp(f1);
    eta2=exp(f2);
    
    nu=1-z;
    sd=lik.stime(2)-lik.stime(1);
    
    if ny==1
      ll=0;
      for i1=1:n
        i3=i3v(i1);
        nft=(i3-1)*ntime;
        si=sum(y(i1)>lik.stime);
        ll=ll + nu(i1).*(f1(si+nft)+f2(i1)) - (y(i1)-lik.stime(si)).*la1(si+nft).*eta2(i1) ...
              - sum(sd.*la1((1+nft):(si-1+nft)).*eta2(i1));
      end
    else
      
      ll=0;
      sb=sum(bsxfun(@gt,y(:,1),lik.stime),2);
      se=sum(bsxfun(@gt,y(:,2),lik.stime),2);
      for i1=1:n        
        i3=i3v(i1);        
        nft=(i3-1)*ntime;
        if sb(i1)==0
          ll=ll + nu(i1).*(f1(se(i1)+nft)+f2(i1)) - (y(i1,2)-lik.stime(se(i1))).*la1(se(i1)+nft).*eta2(i1) ...
                - sum(sd.*la1((1+nft):(se(i1)-1+nft)).*eta2(i1));
        else
          
          if se(i1)==sb(i1)
            ll=ll + nu(i1).*(f1(se(i1)+nft)+f2(i1)) - (y(i1,2)-y(i1,1)).*la1(se(i1)+nft).*eta2(i1);
          else
            ll=ll + nu(i1).*(f1(se(i1)+nft)+f2(i1)) - (y(i1,2)-lik.stime(se(i1))).*la1(se(i1)+nft).*eta2(i1) ...
                - sum(sd.*la1((sb(i1)+1+nft):(se(i1)-1+nft)).*eta2(i1)) - (lik.stime(sb(i1)+1)-y(i1,1)).*la1(sb(i1)+nft).*eta2(i1);
          end
        end
      end
    end
  end

  function llg = lik_coxph_llg(lik, y, f, param, z)
  %LIK_COXPH_LLG  Gradient of the log likelihood
  %
  %  Description 
  %    LLG = LIK_COXPH_LLG(LIK, Y, F, PARAM) takes a likelihood
  %    structure LIK, incedence counts Y, expected counts Z and
  %    latent values F. Returns the gradient of the log likelihood
  %    with respect to PARAM. At the moment PARAM can be 'param' or
  %    'latent'. This subfunction is needed when using Laplace 
  %    approximation or MCMC for inference with non-Gaussian 
  %    likelihoods.
  %
  %  See also
  %    LIK_COXPH_LL, LIK_COXPH_LLG2, LIK_COXPH_LLG3, GPLA_E

    if isempty(z)
      error(['lik_coxph -> lik_coxph_llg: missing z!    '... 
             'Coxph likelihood needs the expected number of    '...
             'occurrences as an extra input z. See, for         '...
             'example, lik_coxph and gpla_e.               ']);
    end
    
    ntime=size(lik.xtime,1);    
    
    [n,ny]=size(y);
    i3v=ones(size(y));
    if isfield(lik, 'stratificationVariables')
      f1=f(1:ntime*lik.n_u);
      f2=f((ntime*lik.n_u+1):(ntime*lik.n_u+n));
      llg=zeros(ntime*lik.n_u+n,1);
      for i=1:lik.n_u
        i3v(lik.stratind{i})=i;
      end
      nf1=lik.n_u*ntime;
    else
      f1=f(1:ntime);
      f2=f((ntime+1):(ntime+n));
      llg=zeros(ntime+n,1);
      nf1=ntime;
    end
    
    la1=exp(f1(:));
    eta2=exp(f2(:));
    
    nu=1-z;
    sd=lik.stime(2)-lik.stime(1);
    
    switch param
      case 'latent'
          
      if ny==1
        if ~isfield(lik, 'stratificationVariables')
          for i1=1:ntime
            ind=y>=lik.stime(i1) & y<lik.stime(i1+1);
            llg(i1)= sum(nu(ind) - (y(ind)-lik.stime(i1)).*la1(i1).*eta2(ind)) - sum(sum(sd.*la1(i1)).*eta2(~ind & y>=lik.stime(i1+1)));
          end
        else
          for j=1:lik.n_u
            for i1=1:ntime
              ind=y>=lik.stime(i1) & y<lik.stime(i1+1);
              ind2=lik.stratind{j};
              llg(i1+(j-1)*ntime)= sum(nu(ind & ind2) - (y(ind & ind2)-lik.stime(i1)).*la1(i1+(j-1)*ntime).*eta2(ind & ind2)) - sum(sum(sd.*la1(i1+(j-1)*ntime)).*eta2((~ind & ind2) & y>=lik.stime(i1+1)));
            end
          end
        end
        for i1=1:n
          i3=i3v(i1);
          nft=(i3-1)*ntime;
          si=sum(y(i1)>lik.stime);
          llg(i1+nf1)= nu(i1) - (y(i1)-lik.stime(si)).*la1(si+nft).*eta2(i1) - sum(sd.*la1((1+nft):(si-1+nft)).*eta2(i1));
        end
      else
        for i1=1:ntime
          
          % left truncated + follow-up entry: (1)
          ind_vkst = y(:,1)>=lik.stime(i1) & y(:,1)<lik.stime(i1+1) & y(:,2)>=lik.stime(i1+1);
          % follow-up entry + follow-up exit: (2)
          ind_stsp = y(:,1)>=lik.stime(i1) & y(:,1)<lik.stime(i1+1) & y(:,2)>=lik.stime(i1) & y(:,2)<lik.stime(i1+1);
          % follow-up: (3)
          ind_s = y(:,1)<lik.stime(i1) & y(:,2)>=lik.stime(i1+1);
          % follow-up exit: (4)
          ind_sp = y(:,1)<lik.stime(i1) & y(:,2)>=lik.stime(i1) & y(:,2)<lik.stime(i1+1);
          
          if ~isfield(lik, 'stratificationVariables')
            % (1)
            s2b=sum(-(lik.stime(i1+1)-y(ind_vkst,1)).*la1(i1).*eta2(ind_vkst));
            % (2)
            s3b=sum(nu(ind_stsp)) - sum((y(ind_stsp,2)-y(ind_stsp,1)).*la1(i1).*eta2(ind_stsp));
            % (3)
            s4= - sum(sd.*la1(i1).*eta2(ind_s));
            % (4)
            s5=sum(nu(ind_sp)) - sum((y(ind_sp,2)-lik.stime(i1)).*la1(i1).*eta2(ind_sp));
            
            llg(i1) = s2b+s3b+s4+s5;
          else
            for j=1:lik.n_u
              ind2 = lik.stratind{j};
              nft=(j-1)*ntime;
              % (1)
              s2b=sum(-(lik.stime(i1+1)-y(ind_vkst & ind2,1)).*la1(i1+nft).*eta2(ind_vkst & ind2));
              % (2)
              s3b=sum(nu(ind_stsp & ind2)) - sum((y(ind_stsp & ind2,2)-y(ind_stsp & ind2,1)).*la1(i1+nft).*eta2(ind_stsp & ind2));
              % (3)
              s4= - sum(sd.*la1(i1+nft).*eta2(ind_s & ind2));
              % (4)
              s5=sum(nu(ind_sp & ind2)) - sum((y(ind_sp & ind2,2)-lik.stime(i1)).*la1(i1+nft).*eta2(ind_sp & ind2));
              
              llg(i1+nft) = s2b+s3b+s4+s5;
            end
          end
          
        end
        
        sb=sum(bsxfun(@gt,y(:,1),lik.stime),2);
        se=sum(bsxfun(@gt,y(:,2),lik.stime),2);
        for i1=1:n
          i3=i3v(i1);
          nft=(i3-1)*ntime;
          if sb(i1)==0
            llg(i1+nf1)= nu(i1) - (y(i1,2)-lik.stime(se(i1))).*la1(se(i1)+nft).*eta2(i1) - sum(sd.*la1((1+nft):(se(i1)-1+nft)).*eta2(i1));
          else
            if se(i1)==sb(i1)
              llg(i1+nf1) = nu(i1) - (y(i1,2)-y(i1,1)).*la1(se(i1)+nft).*eta2(i1);
            else
              llg(i1+nf1) = nu(i1) - (y(i1,2)-lik.stime(se(i1))).*la1(se(i1)+nft).*eta2(i1) - sum(sd.*la1((sb(i1)+1+nft):(se(i1)-1+nft)).*eta2(i1)) - (lik.stime(sb(i1)+1)-y(i1,1)).*la1(sb(i1)+nft).*eta2(i1);
            end
          end
        end
        
      end
      
    end
  end

  function [llg2,llg2mat] = lik_coxph_llg2(lik, y, f, param, z)
  %function [pi_vec, pi_mat] = lik_coxph_llg2(lik, y, ff, param, z)
  %LIK_COXPH_LLG2  Second gradients of the log likelihood
  %
  %  Description        
  %    LLG2 = LIK_COXPH_LLG2(LIK, Y, F, PARAM) takes a likelihood
  %    structure LIK, incedence counts Y, expected counts Z, and
  %    latent values F. Returns the Hessian of the log likelihood
  %    with respect to PARAM. At the moment PARAM can be only
  %    'latent'. LLG2 is a vector with diagonal elements of the
  %    Hessian matrix (off diagonals are zero). This subfunction
  %    is needed when using Laplace approximation or EP for inference 
  %    with non-Gaussian likelihoods.
  %
  %  See also
  %    LIK_COXPH_LL, LIK_COXPH_LLG, LIK_COXPH_LLG3, GPLA_E

    if isempty(z)
      error(['lik_coxph -> lik_coxph_llg2: missing z!   '... 
             'Coxph likelihood needs the expected number of    '...
             'occurrences as an extra input z. See, for         '...
             'example, lik_coxph and gpla_e.               ']);
    end
    
    [n,ny]=size(y);
    ntime=size(lik.xtime,1);
    i3v=ones(size(y,1),1);
    if isfield(lik, 'stratificationVariables')
      f1=f(1:ntime*lik.n_u);
      f2=f((ntime*lik.n_u+1):(ntime*lik.n_u+n));
      llg2=zeros(ntime*lik.n_u+n,1);
      llg2mat=zeros(ntime*lik.n_u,n);
      for i=1:lik.n_u
        i3v(lik.stratind{i})=i;
      end
      nf1=ntime*lik.n_u;
    else
      f1=f(1:ntime);
      f2=f((ntime+1):(ntime+n));
      llg2=zeros(ntime+n,1);
      llg2mat=zeros(ntime,n);
      nf1=ntime;
    end
   
    la1=exp(f1);
    eta2=exp(f2);
    
    %nu=1-z;
    sd=lik.stime(2)-lik.stime(1);
    switch param
      case 'latent'
        
        if ny==1
          % 11
          if ~isfield(lik, 'stratificationVariables')
            for i1=1:ntime
              ind=y>=lik.stime(i1) & y<lik.stime(i1+1);
              llg2(i1)= sum(-(y(ind)-lik.stime(i1)).*la1(i1).*eta2(ind)) - sum(sum(sd.*la1(i1)).*eta2(~ind & y>=lik.stime(i1+1)));
              %llg2(i1,i1)= sum(-(y(ind)-lik.stime(i1)).*la1(i1).*eta2(ind)) - sum(sum(sd.*la1(i1)).*eta2(~ind & y>=lik.stime(i1+1)));
            end
          else
            for j=1:lik.n_u
              for i1=1:ntime                
                ind=y>=lik.stime(i1) & y<lik.stime(i1+1);
                ind2=lik.stratind{j};
                llg2(i1+(j-1)*ntime)= sum(-(y(ind & ind2)-lik.stime(i1)).*la1(i1+ntime*(j-1)).*eta2(ind & ind2)) - sum(sum(sd.*la1(i1+ntime*(j-1))).*eta2((~ind & ind2) & y>=lik.stime(i1+1)));
                %llg2(i1,i1)= sum(-(y(ind)-lik.stime(i1)).*la1(i1).*eta2(ind)) - sum(sum(sd.*la1(i1)).*eta2(~ind & y>=lik.stime(i1+1)));
              end
            end
          end
          
          % 22
          for i1=1:n
            si=sum(y(i1)>lik.stime);
            i3=i3v(i1);
            llg2(i1+nf1)= -(y(i1)-lik.stime(si)).*la1(si+(i3-1)*ntime).*eta2(i1) - sum(sd.*la1((1+(i3-1)*ntime):(si-1+(i3-1)*ntime)).*eta2(i1));
            %llg2(i1+ntime,i1+ntime)= -(y(i1)-lik.stime(si)).*la1(si).*eta2(i1) - sum(sd.*la1(1:(si-1)).*eta2(i1));
          end
          
          % derivative wrt f1 and f2:
          if ~isfield(lik, 'stratificationVariables')
            for i1=1:ntime
              ind=y>=lik.stime(i1) & y<lik.stime(i1+1);
              llg2mat(i1,find(ind))= -(y(ind)-lik.stime(i1)).*la1(i1).*eta2(ind);
              llg2mat(i1,find((~ind & y>=lik.stime(i1+1)))) = - sd.*la1(i1).*eta2((~ind & y>=lik.stime(i1+1)));
              %llg2(i1,ntime+find(ind))= -(y(ind)-lik.stime(i1)).*la1(i1).*eta2(ind);
              %llg2(i1,ntime+find((~ind & y>=lik.stime(i1+1)))) = - sd.*la1(i1).*eta2((~ind & y>=lik.stime(i1+1)));
              %llg2(ntime+find(ind),i1)=llg2(i1,ntime+find(ind));
              %llg2(ntime+find((~ind & y>=lik.stime(i1+1))),i1)=llg2(i1,ntime+find((~ind & y>=lik.stime(i1+1))));
            end
          else
            for j=1:lik.n_u
              ind2=lik.stratind{j};
              for i1=1:ntime
                ind=y>=lik.stime(i1) & y<lik.stime(i1+1);
                llg2mat(i1+(j-1)*ntime,find(ind & ind2))= -(y(ind & ind2)-lik.stime(i1)).*la1(i1+(j-1)*ntime).*eta2(ind & ind2);
                llg2mat(i1+(j-1)*ntime,find(((~ind & ind2) & y>=lik.stime(i1+1)))) = - sd.*la1(i1+(j-1)*ntime).*eta2(((~ind & ind2) & y>=lik.stime(i1+1)));
                %llg2(i1,ntime+find(ind))= -(y(ind)-lik.stime(i1)).*la1(i1).*eta2(ind);
                %llg2(i1,ntime+find((~ind & y>=lik.stime(i1+1)))) = - sd.*la1(i1).*eta2((~ind & y>=lik.stime(i1+1)));
                %llg2(ntime+find(ind),i1)=llg2(i1,ntime+find(ind));
                %llg2(ntime+find((~ind & y>=lik.stime(i1+1))),i1)=llg2(i1,ntime+find((~ind & y>=lik.stime(i1+1))));
              end
            end
          end
          
        else
          
          % 11
          for i1=1:ntime
            
            % left truncated + follow-up entry: (1)
            ind_vkst = y(:,1)>=lik.stime(i1) & y(:,1)<lik.stime(i1+1) & y(:,2)>=lik.stime(i1+1);
            % follow-up entry + follow-up exit: (2)
            ind_stsp = y(:,1)>=lik.stime(i1) & y(:,1)<lik.stime(i1+1) & y(:,2)>=lik.stime(i1) & y(:,2)<lik.stime(i1+1);
            % follow-up: (3)
            ind_s = y(:,1)<lik.stime(i1) & y(:,2)>=lik.stime(i1+1);
            % follow-up exit: (4)
            ind_sp = y(:,1)<lik.stime(i1) & y(:,2)>=lik.stime(i1) & y(:,2)<lik.stime(i1+1);
            
            if ~isfield(lik, 'stratificationVariables')
              % (1)
              s2b=sum(-(lik.stime(i1+1)-y(ind_vkst,1)).*la1(i1).*eta2(ind_vkst));
              % (2)
              s3b=-sum((y(ind_stsp,2)-y(ind_stsp,1)).*la1(i1).*eta2(ind_stsp));
              % (3)
              s4=-sum(sd.*la1(i1).*eta2(ind_s));
              % (4)
              s5=-sum((y(ind_sp,2)-lik.stime(i1)).*la1(i1).*eta2(ind_sp));
              llg2(i1) = s2b+s3b+s4+s5;
            else
              for j=1:lik.n_u
                ind2=lik.stratind{j};
                nft=(j-1)*ntime;
                % (1)
                s2b=sum(-(lik.stime(i1+1)-y(ind_vkst & ind2,1)).*la1(i1+nft).*eta2(ind_vkst & ind2));
                % (2)
                s3b=-sum((y(ind_stsp & ind2,2)-y(ind_stsp & ind2,1)).*la1(i1+nft).*eta2(ind_stsp & ind2));
                % (3)
                s4=-sum(sd.*la1(i1+nft).*eta2(ind_s & ind2));
                % (4)
                s5=-sum((y(ind_sp & ind2,2)-lik.stime(i1)).*la1(i1+nft).*eta2(ind_sp & ind2));
                llg2(i1+nft) = s2b+s3b+s4+s5;
              end
            end
          end
          
          % 22
          sb=sum(bsxfun(@gt,y(:,1),lik.stime),2);
          se=sum(bsxfun(@gt,y(:,2),lik.stime),2);
          for i1=1:n
            i3=i3v(i1);
            nft=(i3-1)*ntime;
            if sb(i1)==0
              llg2(i1+nf1)= -(y(i1,2)-lik.stime(se(i1))).*la1(se(i1)+nft).*eta2(i1) - sum(sd.*la1((1+nft):(se(i1)-1+nft)).*eta2(i1));
              %llg2(i1+ntime,i1+ntime)= -(y(i1,2)-lik.stime(se)).*la1(se).*eta2(i1) - sum(sd.*la1(1:(se-1)).*eta2(i1));
            else
              if se(i1)==sb(i1)
                llg2(i1+nf1) = -(y(i1,2)-y(i1,1)).*la1(se(i1)+nft).*eta2(i1);
                %llg2(i1+ntime,i1+ntime) = -(y(i1,2)-y(i1,1)).*la1(se).*eta2(i1);
              else
                llg2(i1+nf1) = -(y(i1,2)-lik.stime(se(i1))).*la1(se(i1)+nft).*eta2(i1) - sum(sd.*la1((sb(i1)+1+nft):(se(i1)-1+nft)).*eta2(i1)) - (lik.stime(sb(i1)+1)-y(i1,1)).*la1(sb(i1)+nft).*eta2(i1);
                %llg2(i1+ntime,i1+ntime) = -(y(i1,2)-lik.stime(se)).*la1(se).*eta2(i1) - sum(sd.*la1((sb+1):(se-1)).*eta2(i1)) - (lik.stime(sb+1)-y(i1,1)).*la1(sb).*eta2(i1);
              end
            end
          end
          
          % derivative wrt f1 and f2:
          for i1=1:ntime
            
            % left truncated + follow-up entry: (1)
            ind_vkst = y(:,1)>=lik.stime(i1) & y(:,1)<lik.stime(i1+1) & y(:,2)>=lik.stime(i1+1);
            % follow-up entry + follow-up exit: (2)
            ind_stsp = y(:,1)>=lik.stime(i1) & y(:,1)<lik.stime(i1+1) & y(:,2)>=lik.stime(i1) & y(:,2)<lik.stime(i1+1);
            % follow-up: (3)
            ind_s = y(:,1)<lik.stime(i1) & y(:,2)>=lik.stime(i1+1);
            % follow-up exit: (4)
            ind_sp = y(:,1)<lik.stime(i1) & y(:,2)>=lik.stime(i1) & y(:,2)<lik.stime(i1+1);

            
            if ~isfield(lik, 'stratificationVariables')
              % (1)
              llg2mat(i1,ind_vkst)=-(lik.stime(i1+1)-y(ind_vkst,1)).*la1(i1).*eta2(ind_vkst);
              % (2)
              llg2mat(i1,ind_stsp)=-(y(ind_stsp,2)-y(ind_stsp,1)).*la1(i1).*eta2(ind_stsp);
              % (3)
              llg2mat(i1,ind_s)= -sd.*la1(i1).*eta2(ind_s);
              % (4)
              llg2mat(i1,ind_sp)=-(y(ind_sp,2)-lik.stime(i1)).*la1(i1).*eta2(ind_sp);
            else
              for j=1:lik.n_u
                ind2=lik.stratind{j};
                nft=(j-1)*ntime;
                % (1)
                llg2mat(i1+nft,ind_vkst & ind2)=-(lik.stime(i1+1)-y(ind_vkst & ind2,1)).*la1(i1+nft).*eta2(ind_vkst & ind2);
                % (2)
                llg2mat(i1+nft,ind_stsp & ind2)=-(y(ind_stsp & ind2,2)-y(ind_stsp & ind2,1)).*la1(i1+nft).*eta2(ind_stsp & ind2);
                % (3)
                llg2mat(i1+nft,ind_s & ind2)= -sd.*la1(i1+nft).*eta2(ind_s & ind2);
                % (4)
                llg2mat(i1+nft,ind_sp & ind2)=-(y(ind_sp & ind2,2)-lik.stime(i1)).*la1(i1+nft).*eta2(ind_sp & ind2);
              end
            end
          end
      end
    end
  end    
  
  function [llg3,llg3mat] = lik_coxph_llg3(lik, y, f, param, z, j1)
  %LIK_COXPH_LLG3  Third gradients of the log likelihood
  %
  %  Description
  %    LLG3 = LIK_COXPH_LLG3(LIK, Y, F, PARAM) takes a likelihood
  %    structure LIK, incedence counts Y, expected counts Z and
  %    latent values F and returns the third gradients of the log
  %    likelihood with respect to PARAM. At the moment PARAM can be
  %    only 'latent'. LLG3 is a vector with third gradients. This 
  %    subfunction is needed when using Laplace approximation for 
  %    inference with non-Gaussian likelihoods.
  %
  %  See also
  %    LIK_COXPH_LL, LIK_COXPH_LLG, LIK_COXPH_LLG2, GPLA_E, GPLA_G

    if isempty(z)
      error(['lik_coxph -> lik_coxph_llg3: missing z!   '... 
             'Coxph likelihood needs the expected number of    '...
             'occurrences as an extra input z. See, for         '...
             'example, lik_coxph and gpla_e.               ']);
    end
    
    ntime=size(lik.xtime,1);
    
    [n,ny]=size(y);
    if isfield(lik, 'stratificationVariables')
      f1=f(1:ntime*lik.n_u);
      f2=f((ntime*lik.n_u+1):(ntime*lik.n_u+n));
      llg3=zeros(ntime*lik.n_u+n,1);
      llg3mat=zeros(ntime*lik.n_u,n);
      nf1=ntime*lik.n_u;
      if j1>nf1
        for j=1:lik.n_u
          if lik.stratind{j}(j1-nf1)==1;
            i3=j;
            break;
          end
        end
      end
    else
      f1=f(1:ntime);
      f2=f((ntime+1):(ntime+n));
      llg3=zeros(ntime+n,1);
      llg3mat=zeros(ntime,n);
      nf1=ntime;
      i3=1;
    end
    
    la1=exp(f1);
    eta2=exp(f2);
    
    %nu=1-z;
    sd=lik.stime(2)-lik.stime(1);
   
    
    switch param
      case 'latent'
        
        if ny==1
          
          if j1<=nf1
            
            indt=rem(j1,ntime);
            if indt==0
              indt=ntime;
            end
            
            % 11
            ind=y>=lik.stime(indt) & y<lik.stime(indt+1);
            
            if isfield(lik, 'stratificationVariables')
              ind2=lik.stratind{(j1-indt)/ntime+1};
            else
              ind2=ones(size(ind));
            end
            
            fi=find(ind & ind2);
            fni=find((~ind & ind2) & y>=lik.stime(indt));
            
            llg3(j1) = sum(-(y(ind & ind2)-lik.stime(indt)).*la1(j1).*eta2(ind & ind2)) - sum(sum(sd.*la1(j1)).*eta2((~ind & ind2) & y>=lik.stime(indt+1)));

            if ~isempty(fi)
              valtmp=(-(y(ind & ind2)-lik.stime(indt)).*la1(j1).*eta2(ind & ind2));
              for m2i=1:length(valtmp)
                llg3( nf1+fi(m2i))  = valtmp(m2i);
              end
            end
            if ~isempty(fni)
              valtmp2=(-sd.*la1(j1).*eta2(((~ind & ind2) & y>=lik.stime(indt+1))));
              for m2i=1:length(valtmp2)
                llg3( nf1+fni(m2i))  = valtmp2(m2i);
              end
            end
            
            % 12/21
            % derivative wrt f1 and f2:
            val1tmp=-(y(ind & ind2)-lik.stime(indt)).*la1(j1).*eta2(ind & ind2);
            llg3mat(j1,fi)= val1tmp;
            
            val2tmp = - sd.*la1(j1).*eta2(((~ind & ind2) & y>=lik.stime(indt+1)));
            llg3mat(j1,fni) = val2tmp;
            
          else            
            
            % 11
            s1=sum(y(j1-nf1)>lik.stime);
            llg3((1+(i3-1)*ntime):(s1-1+(i3-1)*ntime)) = - sd.*la1((1+(i3-1)*ntime):(s1-1+(i3-1)*ntime)).*eta2(j1-nf1);
            llg3(s1+(i3-1)*ntime) = -(y(j1-nf1)-lik.stime(s1)).*la1(s1+(i3-1)*ntime).*eta2(j1-nf1);
            
            % 22
            llg3(j1) = -(y(j1-nf1)-lik.stime(s1)).*la1(s1+(i3-1)*ntime).*eta2(j1-nf1) ...
                        - sum(sd.*la1((1+(i3-1)*ntime):(s1-1+(i3-1)*ntime)).*eta2(j1-nf1));
            
            % 12/21
            % derivative wrt f1 and f2:
            val3tmp = - sd.*la1((1+(i3-1)*ntime):(s1-1+(i3-1)*ntime)).*eta2(j1-nf1);
            llg3mat((1+(i3-1)*ntime):(s1-1+(i3-1)*ntime),j1-nf1)= val3tmp;
            llg3mat(s1+(i3-1)*ntime,j1-nf1) = -(y(j1-nf1)-lik.stime(s1)).*la1(s1+(i3-1)*ntime).*eta2(j1-nf1);
            
          end
          
        else
          
          if j1<=nf1
            
            indt=rem(j1,ntime);
            if indt==0
              indt=ntime;
            end                      
            if isfield(lik, 'stratificationVariables')
              ind2=lik.stratind{(j1-indt)/ntime+1};
            else
              ind2=ones(size(y,1),1);
            end
            
            % 11
            % left truncated + follow-up entry: (1)
            ind_vkst = y(:,1)>=lik.stime(indt) & y(:,1)<lik.stime(indt+1) & y(:,2)>=lik.stime(indt+1);
            % follow-up entry + follow-up exit: (2)
            ind_stsp = y(:,1)>=lik.stime(indt) & y(:,1)<lik.stime(indt+1) & y(:,2)>=lik.stime(indt) & y(:,2)<lik.stime(indt+1);
            % follow-up: (3)
            ind_s = y(:,1)<lik.stime(indt) & y(:,2)>=lik.stime(indt+1);
            % follow-up exit: (4)
            ind_sp = y(:,1)<lik.stime(indt) & y(:,2)>=lik.stime(indt) & y(:,2)<lik.stime(indt+1);
            
            
            % (1)
            s2b=sum(-(lik.stime(indt+1)-y(ind_vkst & ind2,1)).*la1(j1).*eta2(ind_vkst & ind2));
            % (2)
            s3b=-sum((y(ind_stsp & ind2,2)-y(ind_stsp & ind2,1)).*la1(j1).*eta2(ind_stsp & ind2));
            % (3)
            s4=-sum(sd.*la1(j1).*eta2(ind_s & ind2));
            % (4)
            s5=-sum((y(ind_sp & ind2,2)-lik.stime(indt)).*la1(j1).*eta2(ind_sp & ind2));
            
            llg3(j1) = s2b+s3b+s4+s5;
            
            % 22
            % (1)
            llg3(nf1+find(ind_vkst & ind2))=-(lik.stime(indt+1)-y(ind_vkst & ind2,1)).*la1(j1).*eta2(ind_vkst & ind2);
            % (2)
            llg3(nf1+find(ind_stsp & ind2))=-(y(ind_stsp & ind2,2)-y(ind_stsp & ind2,1)).*la1(j1).*eta2(ind_stsp & ind2);
            % (3)
            llg3(nf1+find(ind_s & ind2))= -sd.*la1(j1).*eta2(ind_s & ind2);
            % (4)
            llg3(nf1+find(ind_sp & ind2))=-(y(ind_sp & ind2,2)-lik.stime(indt)).*la1(j1).*eta2(ind_sp & ind2);
            
            % 12/21
            llg3mat(j1,find(ind_vkst & ind2))=-(lik.stime(indt+1)-y(ind_vkst & ind2,1)).*la1(j1).*eta2(ind_vkst & ind2);
            % (2)
            llg3mat(j1,find(ind_stsp & ind2))=-(y(ind_stsp & ind2,2)-y(ind_stsp & ind2,1)).*la1(j1).*eta2(ind_stsp & ind2);
            % (3)
            llg3mat(j1,find(ind_s & ind2))= -sd.*la1(j1).*eta2(ind_s & ind2);
            % (4)
            llg3mat(j1,find(ind_sp & ind2))=-(y(ind_sp & ind2,2)-lik.stime(indt)).*la1(j1).*eta2(ind_sp & ind2);
          else
            
            sb=sum(y(j1-nf1,1)>lik.stime); % begin
            se=sum(y(j1-nf1,2)>lik.stime); % end
            
            nft=(i3-1)*ntime;
            % 11
            if sb==0
              llg3((1+nft):(se-1+nft))= -sd.*la1((1+nft):(se-1+nft)).*eta2(j1-nf1);
              llg3(se+nft)= -(y(j1-nf1,2)-lik.stime(se)).*la1(se+nft).*eta2(j1-nf1);
            else
              if se==sb
                llg3(se+nft) = -(y(j1-nf1,2)-y(j1-nf1,1)).*la1(se+nft).*eta2(j1-nf1);
              else
                llg3(sb+nft) = - (lik.stime(sb+1)-y(j1-nf1,1)).*la1(sb+nft).*eta2(j1-nf1);
                llg3((sb+1+nft):(se-1+nft)) = - sd.*la1((sb+1+nft):(se-1+nft)).*eta2(j1-nf1);
                llg3(se+nft) = -(y(j1-nf1,2)-lik.stime(se)).*la1(se+nft).*eta2(j1-nf1);
              end
            end
            
            % 12/21
            if sb==0
              llg3mat((1+nft):(se-1+nft),j1-nf1) = -sd.*la1((1+nft):(se-1+nft)).*eta2(j1-nf1);              
              llg3mat(se+nft,j1-nf1)= -(y(j1-nf1,2)-lik.stime(se)).*la1(se+nft).*eta2(j1-nf1);
            else
              if se==sb
                llg3mat(se+nft,j1-nf1) = -(y(j1-nf1,2)-y(j1-nf1,1)).*la1(se+nft).*eta2(j1-nf1);
              else
                llg3mat(sb+nft,j1-nf1) = - (lik.stime(sb+1)-y(j1-nf1,1)).*la1(sb+nft).*eta2(j1-nf1);
                llg3mat((sb+1+nft):(se-1+nft),j1-nf1) = - sd.*la1((sb+1+nft):(se-1+nft)).*eta2(j1-nf1);
                llg3mat(se+nft,j1-nf1) = -(y(j1-nf1,2)-lik.stime(se)).*la1(se+nft).*eta2(j1-nf1);
              end
            end
            
            % 22
            if sb==0
              llg3(j1)= -(y(j1-nf1,2)-lik.stime(se)).*la1(se+nft).*eta2(j1-nf1) - sum(sd.*la1((1+nft):(se-1+nft)).*eta2(j1-nf1));
            else
              if se==sb
                llg3(j1) = -(y(j1-nf1,2)-y(j1-nf1,1)).*la1(se+nft).*eta2(j1-nf1);
              else
                llg3(j1) = -(y(j1-nf1,2)-lik.stime(se)).*la1(se+nft).*eta2(j1-nf1) - sum(sd.*la1((sb+1+nft):(se-1+nft)).*eta2(j1-nf1)) - (lik.stime(sb+1)-y(j1-nf1,1)).*la1(sb+nft).*eta2(j1-nf1);
              end
            end
            
          end
        end
    end
  end

  
  function [logM_0, m_1, sigm2hati1] = lik_coxph_tiltedMoments(lik, y, i1, S2_i, M_i, z)
    %LIK_COXPH_TILTEDMOMENTS  Returns the marginal moments for EP algorithm
    %
    %  Description
    %    [M_0, M_1, M2] = LIK_COXPH_TILTEDMOMENTS(LIK, Y, I, S2,
    %    MYY, Z) takes a likelihood structure LIK, incedence counts
    %    Y, expected counts Z, index I and cavity variance S2 and
    %    mean MYY. Returns the zeroth moment M_0, mean M_1 and
    %    variance M_2 of the posterior marginal (see Rasmussen and
    %    Williams (2006): Gaussian processes for Machine Learning,
    %    page 55). This subfunction is needed when using EP for 
    %    inference with non-Gaussian likelihoods.
    %
    %  See also
    %    GPEP_E
    
      [n,ny]=size(y);
    
      % M_i(end);
      % S2_i(end,end);
      fgrid=M_i(end)+sqrt(S2_i(end,end))*[-6 6];
      fg=linspace(fgrid(1),fgrid(2),15);
      ng=length(fg);
      
      if isfield(gp.lik, 'stratificationVariables')
        ntime=size(lik.xtime,1)*gp.lik.n_u;
      else
        ntime=size(lik.xtime,1);
      end
      
      %f11=f(1:ntime);
      %f2=f((ntime+1):(ntime+n));
      %la1=exp(f1);
      %eta2=exp(f2);
      
      nu=1-z;
      sd=lik.stime(2)-lik.stime(1);
      if ny==1
        sb=1;
        se=sum(bsxfun(@gt,y(i1,1),lik.stime),2);
      end
      indf=sb:se;
      sdvec=[ones(se-1,1)*sd; y(i1)-lik.stime(se)];
      nutmp=zeros(se,1);
      nutmp(se)=nu(i1);
      
      for j1=1:ng
        
        % conditional distribution
        myy=M_i(indf)+S2_i(indf,end)*(1./S2_i(end,end))*(fg(j1)-M_i(end));
        myy0=myy;
        Sigm=S2_i(indf,indf)-S2_i(indf,end)*(1./S2_i(end,end))*S2_i(end,indf);
        Sigm0=Sigm;
        
        nu_prior=Sigm\myy;
        
        nt=size(myy,1);
        c1=exp(fg(j1));
        
        % site parameters
        tautilde=zeros(nt,1);
        nutilde=zeros(nt,1);
        ztilde=zeros(nt,1);
        
        max_small_ep_iter=50;
        tol=1e-9;
        small_ep_iter=1;
        
        tautilde0=Inf; nutilde0=Inf; ztilde0=Inf;
        
        logZep_tmp=0; logZep=Inf;
        %while small_ep_iter <= max_small_ep_iter && (sum(abs(tautilde0-tautilde)>tol) || sum(abs(nutilde0-nutilde)>tol) || sum(abs(ztilde0-ztilde)>tol))
        while small_ep_iter<=max_small_ep_iter && abs(logZep_tmp-logZep)>tol
          logZep_tmp=logZep;
          
          
          %tautilde0=tautilde; nutilde0=nutilde; ztilde0=ztilde;
          
          for k1=1:nt
            
            tau_i=Sigm(k1,k1)^-1-tautilde(k1);
            nu_i = Sigm(k1,k1)^-1*myy(k1)-nutilde(k1);
            myy_i=nu_i/tau_i;
            sigm2_i=tau_i^-1;
            
            % marginal moments
            [logM0(k1), muhati, sigm2hati] = coxph_tiltedMoments(sigm2_i, myy_i, nutmp(k1), sdvec(k1), c1);
            %[M0, muhati, sigm2hati] = coxph_tiltedMoments(lik, y(i1,:), k1, sigm2_i, myy_i, c1, sd_vec(i1), ztmp);
            
            deltatautilde=sigm2hati^-1-tau_i-tautilde(k1);
            tautilde(k1)=tautilde(k1)+deltatautilde;
            nutilde(k1)=sigm2hati^-1*muhati-nu_i;
            
            apu = deltatautilde/(1+deltatautilde*Sigm(k1,k1));
            Sigm = Sigm - apu*(Sigm(:,k1)*Sigm(:,k1)');
            
            % The below is how Rasmussen and Williams
            % (2006) do the update. The above version is
            % more robust.
            %apu = deltatautilde^-1+Sigm(k1,k1);
            %apu = (Sigm(:,k1)/apu)*Sigm(:,k1)';
            %Sigm = Sigm - apu;
            %Sigm=Sigm-(deltatautilde^-1+Sigm(k1,k1))^-1*(Sigm(:,k1)*Sigm(:,k1)');
            
            %myy=Sigm*nutilde;
            myy=Sigm*(nutilde+nu_prior);
            
            muvec_i(k1,1)=myy_i;
            sigm2vec_i(k1,1)=sigm2_i;
            
          end
          
          
          if tautilde > 0
            Stilde=tautilde;
            Stildesqroot=diag(sqrt(tautilde));
            B=eye(nt)+Stildesqroot*Sigm0*Stildesqroot;
            L=chol(B,'lower');
            
            V=(L\Stildesqroot)*Sigm0;
            Sigm=Sigm0-V'*V;
            %myy=Sigm*nutilde;
            myy=Sigm*(nutilde+nu_prior);
            
            %Ls = chol(Sigm);
            
            % Compute the marginal likelihood
            % Direct formula (3.65):
            % Sigmtilde=diag(1./tautilde);
            % mutilde=inv(Stilde)*nutilde;
            %
            % logZep=-0.5*log(det(Sigmtilde+K))-0.5*mutilde'*inv(K+Sigmtilde)*mutilde+
            %         sum(log(normcdf(y.*muvec_i./sqrt(1+sigm2vec_i))))+
            %         0.5*sum(log(sigm2vec_i+1./tautilde))+
            %         sum((muvec_i-mutilde).^2./(2*(sigm2vec_i+1./tautilde)))
            
            % 4. term & 1. term
            term41=0.5*sum(log(1+tautilde.*sigm2vec_i))-sum(log(diag(L)));
            
            % 5. term (1/2 element) & 2. term
            T=1./sigm2vec_i;
            Cnutilde = Sigm0*nutilde;
            L2 = V*nutilde;
            term52 = nutilde'*Cnutilde - L2'*L2 - (nutilde'./(T+Stilde)')*nutilde;
            term52 = term52.*0.5;
            
            % 5. term (2/2 element)
            term5=0.5*muvec_i'.*(T./(Stilde+T))'*(Stilde.*muvec_i-2*nutilde);
            
            % 3. term
            term3 = sum(logM0);
            
            V_tmp=(L\Stildesqroot);
            Sigm_inv_tmp=V_tmp'*V_tmp;
            
            term_add1=-0.5*myy0'*Sigm_inv_tmp*myy0;
            term_add2=myy0'*(eye(nt)-Sigm_inv_tmp*Sigm0)*nutilde;
            logZep = -(term41+term52+term5+term3+term_add1+term_add2);
            
            %logZep = -(term41+term52+term5+term3);
            
            small_ep_iter=small_ep_iter+1;
            %iter=iter+1;
          else
            error('tautilde <= 0')
          end
        end
        
        ZZ(j1,1)=exp(-logZep);
        MM(:,j1)=myy;
        SS2(:,:,j1)=Sigm;
        
      end
      
      %m_0=zeros(1,1);
      %m_1=zeros(nt+1,1);
      %sigm2hati1=zeros(nt+1,nt+1);
      % indf
      
      W=normpdf(fg,M_i(end),sqrt(S2_i(end,end)))*(fg(2)-fg(1));
      
      qw=W.*ZZ';
      m_0=sum(qw);
      m_1=[sum(bsxfun(@times,qw,MM),2); sum(qw.*fg)]./m_0;
      
      m_211=zeros(nt,nt);
      for k1=1:ng
        m_211=m_211+qw(k1)*(SS2(:,:,k1)+MM(:,k1)*MM(:,k1)');
      end
      m_212=(qw.*fg)*MM';
      m_222=(qw.*fg)*fg';
      
      m_2=[m_211 m_212'; m_212 m_222]./m_0;
      
      sigm2hati1=m_2 - m_1*m_1';
      logM_0 = log(m_0);
      
      %figure(1),hold on, plot(fg(j1),logZep,'.')
      %figure(2),hold on, plot(fg(j1),exp(-logZep),'.')
      
  end
  
  function [lpyt,Ey, Vary] = lik_coxph_predy(lik, Ef, Covf, yt, zt)
  %LIK_COXPH_PREDY  Returns the predictive mean, variance and density of y
  %
  %  Description         
  %    [EY, VARY] = LIK_COXPH_PREDY(LIK, EF, VARF) takes a
  %    likelihood structure LIK, posterior mean EF and posterior
  %    Variance VARF of the latent variable and returns the
  %    posterior predictive mean EY and variance VARY of the
  %    observations related to the latent variables. This 
  %    subfunction is needed when computing posterior predictive 
  %    distributions for future observations.
  %        
  %    [Ey, Vary, PY] = LIK_COXPH_PREDY(LIK, EF, VARF YT, ZT)
  %    Returns also the predictive density of YT, that is 
  %        p(yt | zt) = \int p(yt | f, zt) p(f|y) df.
  %    This requires also the incedence counts YT, expected counts ZT.
  %    This subfunction is needed when computing posterior predictive 
  %    distributions for future observations.
  %
  %  See also
  %    GPLA_PRED, GPEP_PRED, GPMC_PRED

    if isempty(zt)
      error(['lik_coxph -> lik_coxph_predy: missing zt!'... 
             'Coxph likelihood needs the expected number of    '...
             'occurrences as an extra input zt. See, for         '...
             'example, lik_coxph and gpla_e.               ']);
    end
    ntime=size(lik.xtime,1);
    
    ntest=size(zt,1);
    ny=size(yt,2);
    
    if isfield(lik, 'stratificationVariables')
      nf1=lik.n_u*ntime;
      i3v=zeros(ntest,1);
      for ii=1:length(lik.stratindt)
        i3v(lik.stratindt{ii})=ii;
      end
    else
      nf1=ntime;
      i3v=ones(ntest,1);
    end
    
    Py = zeros(size(zt));
    %Ey = zeros(size(zt));
    %EVary = zeros(size(zt));
    %VarEy = zeros(size(zt));
    
    S=10000;
    sd=lik.stime(2)-lik.stime(1);
    nu=1-zt;
    
    [nn1,nn2,c] =size(Covf);
    if (c>1) || (nn1==nn2)
      mcmc=false;
    else
      mcmc=true;      
    end
    
    for i1=1:ntest
      i3=i3v(i1);
      if mcmc
        Sigm_tmp=([Covf(1:ntime,:); Covf(i1+ntime,:)]);
        f_star=bsxfun(@plus,Ef([(1+(i3-1)*ntime):(i3*ntime) i1+nf1])', ...
          bsxfun(@times,sqrt(Sigm_tmp'),randn(S,nf1+1)));
      else
        Sigm_tmp=Covf([(1+(i3-1)*ntime):(i3*ntime) i1+nf1],[(1+(i3-1)*ntime):(i3*ntime) i1+nf1]);
        Sigm_tmp=(Sigm_tmp+Sigm_tmp')./2;
        f_star=mvnrnd(Ef([(1+(i3-1)*ntime):(i3*ntime) i1+nf1]), Sigm_tmp, S);
      end
      
      f1=f_star(:,1:ntime);
      f2=f_star(:,(ntime+1):end);
      
      la1=exp(f1);
      eta2=exp(f2);
      
      if ny==1
        si=sum(yt(i1)>lik.stime);
        Py(i1)=mean(exp(nu(i1).*(f1(:,si)+f2) - (yt(i1)-lik.stime(si)).*la1(:,si).*eta2 ...
                     - sum(sd.*la1(:,1:(si-1)),2).*eta2));
      else
        
        sb=sum(bsxfun(@gt,yt(i1,1),lik.stime),2);
        se=sum(bsxfun(@gt,yt(i1,2),lik.stime),2);
        
        if sb==0
          Py(i1) = mean(exp(nu(i1).*(f1(:,se)+f2) - (yt(i1,2)-lik.stime(se)).*la1(:,se).*eta2 - sum(sd.*la1(:,1:(se-1)),2).*eta2));
        else
          if se==sb
            Py(i1) = mean(exp(nu(i1).*(f1(:,se)+f2) - (yt(i1,2)-yt(i1,1)).*la1(:,se).*eta2));
          else
            Py(i1) = mean(exp(nu(i1).*(f1(:,se)+f2) - (yt(i1,2)-lik.stime(se)).*la1(:,se).*eta2 - sum(sd.*la1(:,(sb+1):(se-1)),2).*eta2 - (lik.stime(sb+1)-yt(i1,1)).*la1(:,sb).*eta2));
          end
        end
      end
    end
    Ey = [];
    Vary = [];
    lpyt=log(Py);
    
    %     % Evaluate Ey and Vary
%     for i1=1:length(Ef)
%       %%% With quadrature
%       myy_i = Ef(i1);
%       sigm_i = sqrt(Varf(i1));
%       minf=myy_i-6*sigm_i;
%       maxf=myy_i+6*sigm_i;
% 
%       F = @(f) exp(log(avgE(i1))+f+norm_lpdf(f,myy_i,sigm_i));
%       Ey(i1) = quadgk(F,minf,maxf);
%       
%       F2 = @(f) exp(log(avgE(i1).*exp(f)+((avgE(i1).*exp(f)).^2/r))+norm_lpdf(f,myy_i,sigm_i));
%       EVary(i1) = quadgk(F2,minf,maxf);
%       
%       F3 = @(f) exp(2*log(avgE(i1))+2*f+norm_lpdf(f,myy_i,sigm_i));
%       VarEy(i1) = quadgk(F3,minf,maxf) - Ey(i1).^2;
%     end
%     Vary = EVary + VarEy;
% 
%     % Evaluate the posterior predictive densities of the given observations
%     if nargout > 2
%       for i1=1:length(Ef)
%         % get a function handle of the likelihood times posterior
%         % (likelihood * posterior = Negative-binomial * Gaussian)
%         % and useful integration limits
%         [pdf,minf,maxf]=init_coxph_norm(...
%           yt(i1),Ef(i1),Varf(i1),avgE(i1),r);
%         % integrate over the f to get posterior predictive distribution
%         Py(i1) = quadgk(pdf, minf, maxf);
%       end
%     end
  end

  function [logM_0, m_1, sigm2hati1] = coxph_tiltedMoments(sigm2_i, myy_i, nutmp, sd, c1)
  
  integrand = @(f) exp(-c1.*exp(f).*sd + nutmp*(f+log(c1)) - log(sigm2_i)/2 - log(2*pi)/2 - 0.5*(f-myy_i).^2./sigm2_i);
  RTOL = 1.e-6;
  ATOL = 1.e-10;
  minf=myy_i+sqrt(sigm2_i)*(-6);
  maxf=myy_i+sqrt(sigm2_i)*(6);
  
  [m_0, m_1, m_2] = quad_moments(integrand, minf, maxf, RTOL, ATOL);
  sigm2hati1 = m_2 - m_1.^2;
  
  % If the second central moment is less than cavity variance
  % integrate more precisely. Theoretically for log-concave
  % likelihood should be sigm2hati1 < sigm2_i.
  
  if sigm2hati1 >= sigm2_i
    ATOL = ATOL.^2;
    RTOL = RTOL.^2;
    [m_0, m_1, m_2] = quad_moments(tf, minf, maxf, RTOL, ATOL);
    sigm2hati1 = m_2 - m_1.^2;
    if sigm2hati1 >= sigm2_i
      error('lik_poisson_tilted_moments: sigm2hati1 >= sigm2_i');
    end
  end
  logM_0 = log(m_0);
  end


  function [df,minf,maxf] = init_coxph_norm(yy,myy_i,sigm2_i,avgE,r)
  %INIT_COXPH_NORM
  %
  %  Description
  %    Return function handle to a function evaluating
  %    Negative-Binomial * Gaussian which is used for evaluating
  %    (likelihood * cavity) or (likelihood * posterior) Return
  %    also useful limits for integration. This is private function
  %    for lik_coxph.
  %  
  %  See also
  %    LIK_COXPH_TILTEDMOMENTS, LIK_COXPH_SITEDERIV,
  %    LIK_COXPH_PREDY
    
  % avoid repetitive evaluation of constant part
    ldconst = -gammaln(r)-gammaln(yy+1)+gammaln(r+yy)...
              - log(sigm2_i)/2 - log(2*pi)/2;
    % Create function handle for the function to be integrated
    df = @coxph_norm;
    % use log to avoid underflow, and derivates for faster search
    ld = @log_coxph_norm;
    ldg = @log_coxph_norm_g;
    ldg2 = @log_coxph_norm_g2;

    % Set the limits for integration
    % Negative-binomial likelihood is log-concave so the coxph_norm
    % function is unimodal, which makes things easier
    if yy==0
      % with yy==0, the mode of the likelihood is not defined
      % use the mode of the Gaussian (cavity or posterior) as a first guess
      modef = myy_i;
    else
      % use precision weighted mean of the Gaussian approximation
      % of the Negative-Binomial likelihood and Gaussian
      mu=log(yy/avgE);
      s2=(yy+r)./(yy.*r);
      modef = (myy_i/sigm2_i + mu/s2)/(1/sigm2_i + 1/s2);
    end
    % find the mode of the integrand using Newton iterations
    % few iterations is enough, since the first guess in the right direction
    niter=4;       % number of Newton iterations
    mindelta=1e-6; % tolerance in stopping Newton iterations
    for ni=1:niter
      g=ldg(modef);
      h=ldg2(modef);
      delta=-g/h;
      modef=modef+delta;
      if abs(delta)<mindelta
        break
      end
    end
    % integrand limits based on Gaussian approximation at mode
    modes=sqrt(-1/h);
    minf=modef-8*modes;
    maxf=modef+8*modes;
    modeld=ld(modef);
    iter=0;
    % check that density at end points is low enough
    lddiff=20; % min difference in log-density between mode and end-points
    minld=ld(minf);
    step=1;
    while minld>(modeld-lddiff)
      minf=minf-step*modes;
      minld=ld(minf);
      iter=iter+1;
      step=step*2;
      if iter>100
        error(['lik_coxph -> init_coxph_norm: ' ...
               'integration interval minimun not found ' ...
               'even after looking hard!'])
      end
    end
    maxld=ld(maxf);
    step=1;
    while maxld>(modeld-lddiff)
      maxf=maxf+step*modes;
      maxld=ld(maxf);
      iter=iter+1;
      step=step*2;
      if iter>100
        error(['lik_coxph -> init_coxph_norm: ' ...
               'integration interval maximun not found ' ...
               'even after looking hard!'])
      end
    end
    
%     while minld>(modeld-lddiff)
%       minf=minf-modes;
%       minld=ld(minf);
%       iter=iter+1;
%       if iter>100
%         error(['lik_coxph -> init_coxph_norm: ' ...
%                'integration interval minimun not found ' ...
%                'even after looking hard!'])
%       end
%     end
%     maxld=ld(maxf);
%     while maxld>(modeld-lddiff)
%       maxf=maxf+modes;
%       maxld=ld(maxf);
%       iter=iter+1;
%       if iter>100
%         error(['lik_coxph -> init_coxph_norm: ' ...
%                'integration interval maximum not found ' ...
%                'even after looking hard!'])
%       end
%       
%     end
    
    function integrand = coxph_norm(f)
    % Negative-binomial * Gaussian
      mu = avgE.*exp(f);
      integrand = exp(ldconst ...
                      +yy.*(log(mu)-log(r+mu))+r.*(log(r)-log(r+mu)) ...
                      -0.5*(f-myy_i).^2./sigm2_i);
    end
    
    function log_int = log_coxph_norm(f)
    % log(Negative-binomial * Gaussian)
    % log_coxph_norm is used to avoid underflow when searching
    % integration interval
      mu = avgE.*exp(f);
      log_int = ldconst...
                +yy.*(log(mu)-log(r+mu))+r.*(log(r)-log(r+mu))...
                -0.5*(f-myy_i).^2./sigm2_i;
    end
    
    function g = log_coxph_norm_g(f)
    % d/df log(Negative-binomial * Gaussian)
    % derivative of log_coxph_norm
      mu = avgE.*exp(f);
      g = -(r.*(mu - yy))./(mu.*(mu + r)).*mu ...
          + (myy_i - f)./sigm2_i;
    end
    
    function g2 = log_coxph_norm_g2(f)
    % d^2/df^2 log(Negative-binomial * Gaussian)
    % second derivate of log_coxph_norm
      mu = avgE.*exp(f);
      g2 = -(r*(r + yy))/(mu + r)^2.*mu ...
           -1/sigm2_i;
    end
    
  end

  function p = lik_coxph_invlink(lik, f, z)
  %LIK_COXPH_INVLINK  Returns values of inverse link function
  %             
  %  Description 
  %    P = LIK_COXPH_INVLINK(LIK, F) takes a likelihood structure LIK and
  %    latent values F and returns the values of inverse link function P.
  %    This subfunction is needed when using function gp_predprcty.
  %
  %     See also
  %     LIK_COXPH_LL, LIK_COXPH_PREDY
  
    p = exp(f);
  end



  function [cdf,Ey,Vary] = lik_coxph_predcdf(lik,Ef,Covf,yt)
  %LIK_LOGLOGISTIC_PREDCDF  Returns the predictive cdf evaluated at yt
  %
  %  Description
  %    CDF = LIK_LOGLOGISTIC_PREDCDF(LIK, EF, VARF, YT)
  %    Returns the predictive cdf evaluated at YT given likelihood
  %    structure LIK, posterior mean EF and posterior Variance VARF
  %    of the latent variable. This subfunction is needed when using
  %    functions gp_predcdf or gp_kfcv_cdf.
  %
  %  See also
  %    GP_PREDCDF
  
    
    if isfield(lik, 'stratificationVariables')
      ntime=size(lik.xtime,1)*lik.n_u;
    else
      ntime=size(lik.xtime,1);
    end
    Ef1 = Ef(1:ntime); Ef(1:ntime) = []; Ef2 = Ef;
    nsamps = 10000;
    sd=lik.stime(2)-lik.stime(1);
    Sigm_tmp=Covf;
    [nn1,nn2,cc]=size(Sigm_tmp);
    if cc==1 && nn1~=nn2
      f_star=bsxfun(@plus,[Ef1;Ef2]', ...
        bsxfun(@times,sqrt(Sigm_tmp'),randn(nsamps,ntime+size(Ef2,1))));
    else
      Sigm_tmp=(Sigm_tmp+Sigm_tmp')./2;
      % f_star=mvnrnd(Ef1, Sigm_tmp(1:ntime,1:ntime), nsamps);
      f_star=mvnrnd([Ef1;Ef2], Sigm_tmp, nsamps);
    end
    
    f1=f_star(:,1:ntime);
    f2=f_star(:,(ntime+1):end);
    
    la1=exp(f1);
    eta2=exp(f2);
    
    mST=zeros(size(eta2,2),1);
    if size(yt,2) == 1
      % Integrate from zero to yt
      cumsumtmp=cumsum(la1'*sd)';
      %   t=binsgeq(gp.lik.xtime,yt(i));
      for i1=1:size(eta2,2)
        Stime=exp(-bsxfun(@times,cumsumtmp,eta2(:,i1)));
        mStime=mean(Stime);
        % for i=1:size(yt,1)
        mST(i1)=mStime(binsgeq(lik.xtime,yt(i1)));
        %end
      end
      cdf = 1- mST;
      
    else
      error('Size(yt,2) ~= 1');
    end
    
  end
  
  function reclik = lik_coxph_recappend(reclik, ri, lik)
  %RECAPPEND  Append the parameters to the record
  %
  %  Description 
  %    RECLIK = GPCF_COXPH_RECAPPEND(RECLIK, RI, LIK) takes a
  %    likelihood record structure RECLIK, record index RI and
  %    likelihood structure LIK with the current MCMC samples of
  %    the parameters. Returns RECLIK which contains all the old
  %    samples and the current samples from LIK. This subfunction 
  %    is needed when using MCMC sampling (gp_mc).
  % 
  %  See also
  %    GP_MC

  % Initialize record
    if nargin == 2
      reclik=ri;

      % Set the function handles
      reclik.fh.pak = @lik_coxph_pak;
      reclik.fh.unpak = @lik_coxph_unpak;
      reclik.fh.lp = @lik_coxph_lp;
      reclik.fh.lpg = @lik_coxph_lpg;
      reclik.fh.ll = @lik_coxph_ll;
      reclik.fh.llg = @lik_coxph_llg;    
      reclik.fh.llg2 = @lik_coxph_llg2;
      reclik.fh.llg3 = @lik_coxph_llg3;
      reclik.fh.tiltedMoments = @lik_coxph_tiltedMoments;
      reclik.fh.predy = @lik_coxph_predy;
      reclik.fh.invlink = @lik_coxph_invlink;
      reclik.fh.recappend = @lik_coxph_recappend;
      return
    else
      
      reclik.xtime=lik.xtime;
      reclik.stime=lik.stime;
      
      if isfield(lik, 'stratificationVariables')
        reclik.stratificationVariables=lik.stratificationVariables;
        if isfield(lik,removeStratificationVariables)
          reclik.removeStratificationVariables=lik.removeStratificationVariables;
        end
      end
    end
  end


end

