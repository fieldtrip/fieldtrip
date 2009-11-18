function opt = optimize(obj,data,design,varargin)
%OPTIMIZE optimization function
%
%   This function takes an object, data, and design and optimizes parameters of the object
%   using some validator and evaluation criterion.
%
%    opt = optimize(obj,data,design,varargin)   
%
%   OPTIONS:
%       
%   variables: cell array variable names (fields that occur in obj.options)
%   values: cell array of values that the variables can take on
%   validator: validator object with empty classification procedure
%   criterion: evaluation metric
%   
%   EXAMPLE:
%
%   % optimize function as used by l2svm
%   if isfield(obj.options,'C') && ~isscalar(obj.options.C)              
%     obj.options = optimize(obj,data,design,'variables','C','values',obj.options.C,'validator', ...
%        obj.options.validator,'criterion',obj.options.criterion); 
%   end
%
%  OBSOLETE; use optimizer instead
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: optimize.m,v $
%

    options = parse(varargin{:});

    % the optimizer iterates over all values for the specified
    % variables and returns the best result

    nv = length(options.variables);

    maxresult = -Inf; opt = [];
    
    % iterate over configurations
    configs = cartprod(options.values{:});
    for i=1:size(configs,1)

        % set parameters
        for j=1:nv
            obj.(options.variables{j}) = configs(i,j);
        end

        options.validator.procedure = clfproc({ obj });

        v = options.validator.validate(data,design);

        result = v.evaluate('metric',options.criterion);

        if result > maxresult
            
            maxresult = result;
            opt = obj;
        end

    end
end

function options = parse(varargin)
    
    % parse options
    
    options = struct(varargin{:});

    if ~isfield(options,'validator'), options.validator = crossvalidator(); end
    if ~isfield(options,'verbose'), options.verbose= false; end

    if ~isfield(options,'criterion')
        options.criterion = 'accuracy';
    end

    assert(isa(options.validator,'validator'));

    assert(isfield(options,'variables'));
    assert(isfield(options,'values'));
    if ~iscell(options.variables)
        options.variables = { options.variables };
        options.values = { options.values };
    end

end