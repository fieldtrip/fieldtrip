function CFC =ft_crossfreqanalysis(cfg,LFsig,HFsig)

% input:
% LFsig     low frequency signal to calculate CFC
% HFsig     high frequeny sinal  to calulate CFC
% cfg.method       'plv'-phase locking value;
%                  'mvl'- mean vector length;
%                  'mi' -  modulaiton index 
%cfg.keeptrials    'no' - averaged CFC
%                  'yes'-  calculate CFC each trial

% output:
% cfcs            cross frequency coupling
% LF              low frequency vector
% HF              High frequency vector

LF              = LFsig.freq;
HF              = HFsig.freq;
[ntrial,nchann,~,~] = size(LFsig.fourierspctrm);

% prepare cfcdata 
switch cfg.method
    
    case 'plv'         % phase locking value
       
    plvdatas = zeros(ntrial,nchann,numel(LF),numel(HF)) ;    
    
     for  i =1:nchann       
          chandataLF = LFsig.fourierspctrm(:,i,:,:);
          chandataHF = HFsig.fourierspctrm(:,i,:,:);   
     for j = 1:ntrial        
      plvdatas(j,i,:,:) = data2plv(squeeze(chandataLF(j,:,:,:)),squeeze(chandataHF(j,:,:,:)));       
     end     
     end
     
     cfcdata  = plvdatas;
        
    case  'mvl'   % mean vector length
        
     mvldatas = zeros(ntrial,nchann,numel(LF),numel(HF));     
    
     for  i =1:nchann       
          chandataLF = LFsig.fourierspctrm(:,i,:,:);
          chandataHF = HFsig.fourierspctrm(:,i,:,:);     
     for j = 1:ntrial        
      mvldatas(j,i,:,:) = data2mvl(squeeze(chandataLF(j,:,:,:)),squeeze(chandataHF(j,:,:,:)));       
     end     
     end
     
     cfcdata  = mvldatas;
        
        
    case  'mi'  %  modulaiton index
     
    nbin        = 20;                      % number of phase bin 
    pacdatas   = zeros(ntrial,nchann,numel(LF),numel(HF),nbin) ; 
     
    
     for  i =1:nchann       
          chandataLF = LFsig.fourierspctrm(:,i,:,:);
          chandataHF = HFsig.fourierspctrm(:,i,:,:);  
     for j = 1:ntrial        
      pacdatas(j,i,:,:,:) = data2pac(squeeze(chandataLF(j,:,:,:)),squeeze(chandataHF(j,:,:,:)),nbin);       
     end     
     end
     
     cfcdata  = pacdatas;
    
end
  
 methods.method      = cfg.method ;
 methods.keeptrials  = cfg.keeptrials;
 
 % cfc qualification
 
 [cfcs,dimord] = cfc_quantify(cfcdata,methods);         
        

CFC.cfcs    =cfcs;
CFC.dimord  = dimord;
CFC.LF      = LF;
CFC.HF      = HF; 




end


function [plvdata] =data2plv(LFsigtemp,HFsigtemp)
   
    LFphas   = angle(LFsigtemp);    
    HFamp    = abs(HFsigtemp);
    HFamp(isnan(HFamp(:))) = 0;                  % replace nan with 0
    HFphas   = angle(hilbert(HFamp'))';
    plvdata  = zeros(size(LFsigtemp,1),size(HFsigtemp,1));   % phase lokcing value 

   for i =  1:size(LFsigtemp)
    for j = 1:size(HFsigtemp)
    plvdata(i,j) = nanmean(exp(1i*(LFphas(i,:)-HFphas(j,:)))); 
    end
   end

end

function [mvldata] =data2mvl(LFsigtemp,HFsigtemp)
 % calculate  mean vector length (complex value) per trial 
 % mvldata dim: LF*HF
   
    LFphas   = angle(LFsigtemp);    
    HFamp    = abs(HFsigtemp);   
    mvldata  = zeros(size(LFsigtemp,1),size(HFsigtemp,1));    % mean vector length 
 
   for i =  1:size(LFsigtemp)
   for j = 1:size(HFsigtemp)    
    mvldata(i,j) = nanmean(HFamp(j,:).*exp(1i*LFphas(i,:)));
   end
   end

end


function  pacdata =data2pac(LFsigtemp,HFsigtemp,nbin)
 
 % calculate phase amplitude distribution per trial 
 % pacdata dim: LF*HF*Phasebin
 pacdata  = zeros(size(LFsigtemp,1),size(HFsigtemp,1),nbin);
 
     Ang  = angle(LFsigtemp);
     Amp  = abs(HFsigtemp);
    [~,bin] = histc(Ang, linspace(-pi,pi,nbin));  % binned low frequency phase    
    binamp = zeros (size(HFsigtemp,1),nbin);          %  binned amplitude 
    
    for i= 1:size(Ang,1)
    for k =1:nbin
      idx = find(bin(i,:)==k);
      binamp(:,k)   = mean(Amp(:,idx),2);    
    end
   pacdata(i,:,:) = binamp;    
    end    
end


function [cfcs,dimord] = cfc_quantify(cfcdata,methods)

switch methods.method
    
    case 'plv'  
      
      if strcmp(methods.keeptrials,'no')
      cfcs   = squeeze(abs(mean(cfcdata,1)));
      dimord = 'chan_lf_hf' ;
      else
      cfcs   = abs(cfcdata);
      dimord = 'rpt_chan_lf_hf' ;
      end        
        
    case  'mvl'
        
     if strcmp(methods.keeptrials,'no')
      cfcs   = squeeze(abs(mean(cfcdata,1)));
      dimord = 'chan_lf_hf' ;
      else
      cfcs   = abs(cfcdata);
      dimord = 'rpt_chan_lf_hf' ;
      end     
        
    case  'mi'
        
      [ntrial,nchan,nlf,nhf,nbin] = size(cfcdata);
      
     if strcmp(methods.keeptrials,'yes')
       dimord = 'rpt_chan_lf_hf' ;   
       cfcs = zeros(ntrial,nchan,nlf,nhf);
       for k =1:ntrial
        for n=1:nchan
       pac = squeeze(cfcdata(k,n,:,:,:));      
       Q =ones(nbin,1)/nbin;  % uniform distribution
       mi = zeros(nlf,nhf);
       
  for i=1:nlf
  for j=1:nhf
  P = squeeze(pac(i,j,:))/ nansum(pac(i,j,:));   % normalized distribution 
  % KL distance 
  mi(i,j) = nansum(P.* (log2(P)-log2(Q)))/log2(pha);
  end
  end
  cfcs(k,n,:,:) =mi;  
  
        end
       end    
    
    
     else
     dimord = 'chan_lf_hf' ;   
     cfcs = zeros(nchan,nlf,nhf);   
     cfcdatamean = squeeze(mean(cfcdata,1));
     
     for k =1:nchan
      pac = squeeze(cfcdatamean(k,:,:,:));
      Q =ones(nbin,1)/nbin;                      % uniform distribution
      mi = zeros(nlf,nhf); 
           
  for i=1:nlf
  for j=1:nhf
  P = squeeze(pac(i,j,:))/ nansum(pac(i,j,:));   % normalized distribution 
  % KL distance 
  mi(i,j) = nansum(P.* (log2(P)-log2(Q)))/log2(nbin);
  end
  end
      cfcs(k,:,:) = mi;
     end
     
         
     end
     

end

end
 
 
 