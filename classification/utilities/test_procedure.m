function [acc,sig,cv] = test_procedure(myproc,cvfolds,data,design)
% test classification procedure
%
% [acc,sig,cv] = test_procedure(myproc,cvfolds,data,design)
%
%   Copyright (c) 2009, Marcel van Gerven
%
%   $Log: test_procedure.m,v $
%

  if nargin == 1
    cvfolds = 10;
  end

  if nargin < 4
    
    % load example data
    load covattfrq1
    
    cvec = ismember(left.label,channelselection({'MLO' 'MRO'},left.label)); % subset of channels
    fvec = (left.freq >= 8 & left.freq <= 14); % subset of frequencies
    tvec = (left.time >= 1.5 & left.time <= 2.5); % subset of time segment
    
    data    = [squeeze(mean(left.powspctrm(:,cvec,fvec,tvec),4)); squeeze(mean(right.powspctrm(:,cvec,fvec,tvec),4))];
    design  = [ones(size(left.powspctrm,1),1); 2*ones(size(right.powspctrm,1),1)];
  
    if isa(myproc{end},'transfer_learner')
      
      load covattfrq2;
      
      cvec = ismember(left.label,channelselection({'MLO' 'MRO'},left.label)); % subset of channels
      fvec = (left.freq >= 8 & left.freq <= 14); % subset of frequencies
      tvec = (left.time >= 1.5 & left.time <= 2.5); % subset of time segment
      
      data2    = [squeeze(mean(left.powspctrm(:,cvec,fvec,tvec),4)); squeeze(mean(right.powspctrm(:,cvec,fvec,tvec),4))];
      design2  = [ones(size(left.powspctrm,1),1); 2*ones(size(right.powspctrm,1),1)];
      
      data = {data data2};
      design = {design design2};
      
    end    
    
  end
  
  cv = crossvalidator('procedure',myproc,'cvfolds',cvfolds,'randomize',true,'verbose',true,'compact',false,'balanced',false);

  if isa(myproc{end},'transfer_learner') && ~iscell(data)
    cv = cv.validate({data data},{design design});
  else
    cv = cv.validate(data,design);
  end
  
  acc = cv.evaluate;
  sig = cv.significance;
  
end