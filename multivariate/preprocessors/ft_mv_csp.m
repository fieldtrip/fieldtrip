classdef ft_mv_csp < ft_mv_preprocessor
%FT_MV_CSP extracts and applies CSP (or EED) filters
%
%   Properties - options:
%   
%   'numchan'       : the number of channels (the second data dim);
%                     estimated from data if unspecified
%   'numpatterns'   : number of CSP patterns to extract (2*numpatterns vector filters)
%   'filttype'      : type of filter to apply ('stdCOV','EED','CSP0','CSP3')
%                     'CSP0' ->  CSP by common diagonalisation of class1+class2 
%                                covariance matrix (biosig-> csp.m, see also csp_wrapper.m)
%                     'CSP3' ->  CSP calculation as generalized eigenvalues
%                                (biosig-> csp.m, see also csp_wrapper.m)
%                     'stdCOV'-> based on standard Fieldtrip estimation of
%                                covariance matrix 
%                     'EED'   -> extreme energy difference criterion
%   'outputdatatype': specifies the format of the output data (after testing)
%      'rawcsp' corresponds to data processed using csp filters (time courses) - the output is in the format corresponding to that of the input
%      'powcsp' corresponds to power estimates evaluted as the sum of square of csp-filtered samples - the output is a feature matrix with observations in rows
%      'logpowcsp' is log10(powcsp) - the output is a feature matrix with observations in rows
%
% 
% PARAMETERS:
%  filters;      % extracted filters
%  eigenvalues;  % corresponding eigenvalues (or singular values)
%
%   REMARK
%           Input data should be specified as 3-dimensional array,
%           e.g. [observations(trials) x chan x time].In most applications,
%           this 3-dimensional array can be obtained as a result of TIMELOCKANALYSIS.
%           Example:
%                   cfg = [];
%                   cfg.keeptrials = 'yes'; 
%                   tlcond1 = timelockanalysis(cfg,cond1_data);
%                   tlcond2 = timelockanalysis(cfg,cond2_data);
%                   data    = [tlcond1.trial; tlcond2.trial];
%                   design  = [1*ones(size(tlcond1.trial,1),1); 2*ones(size(tlcond2.trial,1),1)];
%
%
%   EXAMPLE
%    Linear classification with CSP feature extraction 
%
%    myproc = mva({cspprocessor('numchan',151,'numpatterns',3,'outputdatatype','logpowcsp','filttype','CSP0')...
%                      da()});
%
%   SEE ALSO
%   csp_train.m
%   csp_test.m
%
%   REFERENCES
%   Common Spatial Patterns (CSP) - 
%                Ramoser H., Gerking J.M., Pfurtscheller G., Optimal spatial filtering of single trial EEG during 
%                imagined hand movement. IEEE trans Rehab. Eng.,vol. 8, pp.441-446, 2000.    
%                 
%                Blankertz B., Curio G., Muller  K.-R., Classifying Single Trial EEG:
%                Towards Brain Computer Interfacing, in: T. G. Diettrich, S. Becker, and
%                Z. Ghahramani, eds., Advances in Neural Inf. Proc. Systems (NIPS 01),vol. 14, 157�164, 2002.
%
%
%   Extreme Energy Difference (EED) criterion -  
%                Li J.,Sun S., Energy feature extraction of EEG signals and a case study, IJCNN 2008
%                
%                Sun S., The Extreme Energy Ratio Criterion for EEG Feature Extraction, Lecture Notes in Computer Science, vol. 5164/2008
%
%
%   Copyright (c) 2008, Pawel Herman, Marcel van Gerven


    properties        

      numchan
      numpatterns = 1;
      filttype = 'CSP0';
      outputdatatype = 'powcsp';  %'rawcsp' or 'logpowcsp'
       
      filters
      eigenvalues
      
    end

    methods
      
      function obj = ft_mv_csp(varargin)
        
        obj = obj@ft_mv_preprocessor(varargin{:});
        
      end
      
      function obj = train(obj,X,Y)
                
        % multiple datasets
        if iscell(X) || iscell(Y)
          obj = ft_mv_ndata('mvmethod',obj);
          obj = obj.train(X,Y);
          return;
        end
        
        % missing data
        if any(isnan(X(:))) || any(isnan(Y(:))), error('method does not handle missing data'); end
        
        % multiple outputs
        if size(Y,2) > 1
          obj = ft_mv_noutput('mvmethod',obj);
          obj = obj.train(X,Y);
          return;
        end
        
        if isempty(obj.numchan)
          error('number of channels not specified');
        end
        
        if isnumeric(X) && length(size(X)) == 2  %just to make it explicit
          X = reshape(X,size(X,1),obj.numchan,size(X,2)/obj.numchan);
        elseif isnumeric(data) && length(size(X)) == 3   %only for simulating outside mva pipe
          warning('The object used outside any CLFPROC pipe'); %#ok<WNTAG>
        else
          error('Unknown data format for training');
        end
        
        [obj.filters,obj.eigenvalues] = csp_train(X,Y,obj.numpatterns,obj.filttype);
      
      end
      
      function Y = test(obj,X)
        
        if isnumeric(X) && length(size(X)) == 2  %just to make it explicit
          X = reshape(X,size(X,1),obj.numchan,size(X,2)/obj.numchan);
        elseif isnumeric(X) && length(size(X)) == 3   %only for simulating outside mva pipe
          warning('The object used outside any CLFPROC pipe'); %#ok<WNTAG>
        else
          error('Unknown data format for testing');
        end
        
        if strcmp(obj.outputdatatype,'rawcsp')
          Y = csp_test(X,obj.filters);
        elseif strcmp(obj.outputdatatype,'powcsp')
          [aux,csp_pow] = csp_test(X,obj.filters);
          Y = csp_pow;
        elseif strcmp(obj.outputdatatype,'logpowcsp')
          [aux,csp_pow] = csp_test(X,obj.filters);
          Y = log10(csp_pow);
        else
          error('Output data type is not specified - it can be rawcsp,powcsp or logpowcsp');
        end
        
      end
      
    end
end
