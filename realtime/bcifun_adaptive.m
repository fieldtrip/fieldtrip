function cmd = bcifun_adaptive(cfg,data)
% BCIFUN_ADAPTIVE either computes delta rule or error potential output
% depending on triggers. It uses the following cfg. Note that the
% persistent variables need to be cleared with 'clear bcifun_adaptive' 
% whenever the experiment is restarted!
% 
% cfg = 
% 
% triggers: [1 2]
% blocksize: 1
% offset: 0
% bufferdata: 'first'
% bcifun: @bcifun_adaptive
% readevent: 'yes'
% dataset: '/Volumes/Elements/data/alphalat/subject11/Sander_AlphaLat_20080530_01.ds'
% channel: {'MLO'  'MRO'}
%
% for online we need:
%      bufferdata: 'last'
%      dataset: 'shm://'
%      ostream: 'tcp://presentation011:1976'
%
% Copyright (C) 2009, Marcel van Gerven, Ali Bahramisharif, Vicenc Gomez, Alberto LLera

persistent cfgpl; % planar stuff
persistent cfgf1; % stuff for covert attention classifier
persistent cfgf2; % stuff for error potential classifier

sigm = @(x)(1./(1+exp(-x)));
deriv_sigm = @(x)(sigm(x).*(1-sigm(x)));

% loop over acquired trials
cmd = nan(1,length(data.trial));
for trllop=1:length(data.trial)
  
  if isempty(cfgpl)
    
    load grad;
    
    cfgpl.channel      = 'all';
    cfgpl.planarmethod = 'sincos';
    cfgpl.grad = grad;
    cfgpl.online=1;
    
  end
  
  data.grad=cfgpl.grad;
  plan = megplanar(cfgpl, data);
  
  if ~isfield(cfgpl,'onlineprocess')
    cfgpl.onlineprocess=plan.cfg.onlineprocess;
  end
  
  if data.cfg.trl(trllop,4) == 1 % covert attention classifier (trigger = 1)
    
    if isempty(cfgf1)
      
      cfgf1.output       = 'pow';
      cfgf1.channel      = 'all';%{'MLO33_dH','MLO33_dV','MRO33_dH','MRO33_dV'};
      cfgf1.keeptrials   = 'yes';
      cfgf1.method       = 'mtmconvol';
      cfgf1.foi          = 10;
      cfgf1.toi          = 0.5 * cfg.blocksize; % FIXME: Ali
      cfgf1.t_ftimwin    = ones(1,length(cfgf1.foi)) * (cfg.blocksize-0.001); % 0.5 s. timewindow
      cfgf1.taper        = 'hanning';
      cfgf1.online=1;
      
    end
    
    comb = freqanalysis(cfgf1,plan);
    
    if ~isfield(cfgf1,'onlineprocess') || ~isfield(cfgf1.onlineprocess,'cfg')
      
      cfgf1.onlineprocess=comb.cfg.onlineprocess;
      
      %planar gradient
      comb.grad.type='ctf275_planar';
      planar    = planarchannelset(comb);
      cfgf1.sel_dH    = match_str(comb.label, planar(:,1));  % indices of the horizontal channels
      cfgf1.sel_dV    = match_str(comb.label, planar(:,2));  % indices of the vertical  channels
      [dum, sel_planar] = match_str(comb.label, planar(:,1));
      comb.label=planar(sel_planar,3);
      
    end
    
    comb.powspctrm=comb.powspctrm(:,cfgf1.sel_dH,:) + comb.powspctrm(:,cfgf1.sel_dV,:);
    
    if ~isfield(cfgf1,'sel_L')
      %%FIXME: just 'MLO' or 'MRO' is not working!
      listL={'MLO33'};%{'MLO11'    'MLO13'    'MLO14'    'MLO21'    'MLO22'    'MLO23'   'MLO24' 'MLO31'   'MLO32'    'MLO33'   'MLO34'    'MLO41'    'MLO42'    'MLO43'  'MLO44'    'MLO51' 'MLO52'    'MLO53'};
      listR={'MRO33'};%{ 'MRO11' 'MRO12' 'MRO13'    'MRO14' 'MRO21'    'MRO22'    'MRO23'    'MRO24'    'MRO31'    'MRO32'    'MRO33' 'MRO34' 'MRO41'    'MRO42'    'MRO43'    'MRO44'    'MRO51'    'MRO52' 'MRO53'};
      cfgf1.sel_L=match_str(comb.label,channelselection(listL,comb.label));
      cfgf1.sel_R=match_str(comb.label,channelselection(listR,comb.label));
    end
    
    %the parameter of interest is calculated in the next line
    x=log10(mean(comb.powspctrm(cfgf1.sel_L))./mean(comb.powspctrm(cfgf1.sel_R)));
    
    % output of adaptive classifier C1
    if ~isfield(cfgf1,'w0')
      cfgf1.w0 = -1;
      cfgf1.w = 1/100;
    end
    
    cfgf1.prevclf = sigm((x+100)*cfgf1.w + cfgf1.w0);
    cfgf1.prevali = x;
    
    cmd(trllop) = round(cfgf1.prevclf) + 1;
    
    
  else % error potential classifier (trigger = 2)
    
    if isempty(cfgf2)
      
      cfgf2.output       = 'pow';
      cfgf2.channel      = 'all';%{'MLO33_dH','MLO33_dV','MRO33_dH','MRO33_dV'};
      cfgf2.keeptrials   = 'yes';
      cfgf2.method       = 'mtmconvol';
      cfgf2.foi          = 8:30;
      cfgf2.toi          = 0.5 * cfg.blocksize; % FIXME: Ali
      cfgf2.t_ftimwin    = ones(1,length(cfgf2.foi)) * (cfg.blocksize-0.001); % 0.5 s. timewindow
      cfgf2.taper        = 'hanning';
      cfgf2.online=1;
      
    end
    
    comb = freqanalysis(cfgf2,plan);
    
    if ~isfield(cfgf2,'onlineprocess') || ~isfield(cfgf2.onlineprocess,'cfg')
      
      cfgf2.onlineprocess=comb.cfg.onlineprocess;
      
      %planar gradient
      comb.grad.type='ctf275_planar';
      planar    = planarchannelset(comb);
      cfgf2.sel_dH    = match_str(comb.label, planar(:,1));  % indices of the horizontal channels
      cfgf2.sel_dV    = match_str(comb.label, planar(:,2));  % indices of the vertical  channels
      [dum, sel_planar] = match_str(comb.label, planar(:,1));
      comb.label=planar(sel_planar,3);
      
    end
    
    comb.powspctrm=comb.powspctrm(:,cfgf2.sel_dH,:) + comb.powspctrm(:,cfgf2.sel_dV,:);
    
    if ~isfield(cfgf2,'sv')
      
      load params;
      
      cfgf2.sv = sv;
      cfgf2.weights = weights;
      cfgf2.b = b;
      cfgf2.kernel = kernel;
      cfgf2.kerparam = kerparam;

    end
    
    % FIXME: Alberto's responsibility
    cmd(trllop) = svm_eval(comb.powspctrm,cfgf2.sv,cfgf2.weights,cfgf2.b,cfgf2.kernel,cfgf2.kerparam);
    
    if cmd(trllop) == 1 % update delta rule in case of an error
      
      cfgf1.w = cfgf1.w + 1e-5 * ((round(1-cfgf1.prevclf)+1) - cfgf1.prevclf) * deriv_sigm(cfgf1.prevclf) * cfgf1.prevali;
      
    end
    
  end
  
  first = false;
  
end % loop over trials