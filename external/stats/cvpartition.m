%CVPARTITION Create a cross-validation partition for data.
%   An object of the CVPARTITION class defines a random partition on a
%   set of data of a specified size.  This partition can be used to
%   define test and training sets for validating a statistical model
%   using cross-validation.
%
%   C = CVPARTITION(N,'KFold',K) creates a CVPARTITION object C
%   defining a random partition for K-fold cross-validation on N
%   observations. The partition divides N observations into K disjoint
%   subsamples (folds), chosen randomly but with roughly equal size.
%   The default value of K is 10.
%
%   C = CVPARTITION(GROUP,'KFold',K) creates a CVPARTITION object C
%   defining a random partition for a stratified K-fold
%   cross-validation. GROUP is a vector indicating the class
%   information for each observation. GROUP can be a categorical
%   variable, a numeric vector, a string array, or a cell array of
%   strings. Each subsample has roughly equal size and roughly the same
%   class proportions as in GROUP. CVPARTITION treats NaNs or empty
%   strings in GROUP as missing values.
%
%   C = CVPARTITION(GROUP,'KFold',K,'Stratify',stratifyOption) returns
%   an object C defining a random partition for k-fold cross
%   validation. If you supply group as the first input to cvpartition,
%   the function implements stratification by default. If,
%   additionally, you specify 'Stratify',false, the function ignores
%   the class information in group and creates nonstratified random
%   partitions. You cannot specify 'Stratify',true when the first input
%   to cvpartition is n.
%
%   C = CVPARTITION(N,'HoldOut',P) creates a CVPARTITION object C
%   defining a random partition for hold-out validation on N
%   observations. This partition divides N observations into a training
%   set and a test set. P must be a scalar. When 0<P<1, CVPARTITION
%   randomly selects approximately P*N observations for the test set.
%   When P is an integer, it randomly selects P observations for the
%   test set. The default value of P is 1/10.
%
%   C = CVPARTITION(GROUP,'HoldOut',P) randomly partitions observations
%   into a training set and a test set with stratification, using the
%   class information in GROUP, i.e., both training and test sets have
%   roughly the same class proportions as in GROUP.
%
%   C = CVPARTITION(GROUP,'HoldOut',P,'Stratify',stratifyOption)
%   returns an object C defining a random partition into a training set
%   and a holdout or test set. If you supply group as the first input
%   to cvpartition, the function implements stratification by default.
%   If, additionally, you specify 'Stratify',false, the function
%   ignores the class information in group and creates nonstratified
%   random partitions. You cannot specify 'Stratify',true when the
%   first input to cvpartition is n.
%
%   C = CVPARTITION(N,'LeaveOut') creates an object C defining a random
%   partition for leave-one-out cross-validation on N observations.
%   Leave-one-out is a special case of K-fold in which the number of
%   folds is equal to the number of observations N.
%
%   C = CVPARTITION(N,'Resubstitution') creates a CVPARTITION object C
%   which does not partition the data. Both the training set and the
%   test set contain all of the original N observations.
%
%   C = CVPARTITION('CustomPartition',TestSets) creates a CVPARTITION
%   object C which partitions the data based on the test sets indicated
%   in TestSets.
%
%   CVPARTITION properties:
%      Type               - Type of validation partition.
%      NumObservations    - Number of observations.
%      NumTestSets        - Number of test sets.
%      TrainSize          - Size of each training set.
%      TestSize           - Size of each test set.
%      IsCustom           - Indicator for custom partition.
%
%   CVPARTITION methods:
%      repartition        - Rerandomize a cross-validation partition.
%      test               - Test set for a cross-validation partition.
%      training           - Training set for a cross-validation partition.
%
%   Example: Use 10-fold stratified cross-validation to compute the
%   misclassification error for CLASSIFY on iris data.
%
%     load('fisheriris');
%     CVO = cvpartition(species,'KFold',10);
%     err = zeros(CVO.NumTestSets,1);
%     for i = 1:CVO.NumTestSets
%         trIdx = CVO.training(i);
%         teIdx = CVO.test(i);
%         ytest = classify(meas(teIdx,:),meas(trIdx,:),species(trIdx,:));
%         err(i) = sum(~strcmp(ytest,species(teIdx)));
%     end
%     cvErr = sum(err)/sum(CVO.TestSize);
%
%   See also CVPARTITION/TEST, CVPARTITION/TRAINING, CVPARTITION/REPARTITION, CROSSVAL.

%   Copyright 2026 Jan-Mathijs Schoffelen, with the help of Mistral.ai
%     The docstring is from the original matlab version, but the code is a rewrite

classdef cvpartition
  properties
    Type
    NumObservations
    NumTestSets
    TrainSize
    TestSize
    IsCustom
    TestSets
    Group
  end

  methods
    function obj = cvpartition(varargin)
      % Constructor for cvpartition
      if nargin == 0
        error('Not enough input arguments.');
      end

      % Parse input arguments
      if ischar(varargin{1}) && strcmpi(varargin{1}, 'CustomPartition')
        obj = obj.createCustomPartition(varargin{2});
      elseif isnumeric(varargin{1}) && isscalar(varargin{1})
        % First input argument is a scalar
        obj = obj.createNumericPartition(varargin{:});
      elseif iscell(varargin{1}) || isstring(varargin{1}) || iscategorical(varargin{1}) || isnumeric(varargin{1})
        % First input argument is an array that indexes the class for each of the observations
        varargin{1} = varargin{1}(:); % ensure column
        obj = obj.createGroupPartition(varargin{:});
      else
        error('Invalid input arguments, the first argument should either be a scalar, an array that indexes group-membership of the observations, or ''CustomPartition''.');
      end
    end

    function obj = createNumericPartition(obj, N, varargin)
      obj.NumObservations = N;
      obj.IsCustom        = false;

      % Default values
      partitionType = 'KFold';
      K = 10;
      P = 0.1;

      % Parse varargin
      i = 1;
      while i <= length(varargin)
        arg = varargin{i};
        if ischar(arg)
          switch lower(arg)
            case 'kfold'
              partitionType = 'KFold';
              K = varargin{i+1};
              i = i + 2;
            case 'holdout'
              partitionType = 'HoldOut';
              P = varargin{i+1};
              i = i + 2;
            case 'leaveout'
              partitionType = 'LeaveOut';
              i = i + 1;
            case 'resubstitution'
              partitionType = 'Resubstitution';
              i = i + 1;
            otherwise
              error('Unknown partition type: %s', arg);
          end
        else
          error('Invalid argument type.');
        end
      end

      % Create partition
      switch partitionType
        case 'KFold'
          obj = obj.createKFold(N, K);
        case 'HoldOut'
          obj = obj.createHoldOut(N, P);
        case 'LeaveOut'
          obj = obj.createLeaveOut(N);
        case 'Resubstitution'
          obj = obj.createResubstitution(N);
      end
    end

    function obj = createGroupPartition(obj, group, varargin)
      obj.Group           = group;
      obj.NumObservations = length(group);
      obj.IsCustom        = false;

      % Default values
      partitionType = 'KFold';
      K = 10;
      P = 0.1;
      stratifyOption = true;

      % Parse varargin
      i = 1;
      while i <= length(varargin)
        arg = varargin{i};
        if ischar(arg)
          switch lower(arg)
            case 'kfold'
              partitionType = 'KFold';
              K = varargin{i+1};
              i = i + 2;
            case 'holdout'
              partitionType = 'HoldOut';
              P = varargin{i+1};
              i = i + 2;
            case 'stratify'
              stratifyOption = varargin{i+1};
              i = i + 2;
            otherwise
              error('Unknown partition type: %s', arg);
          end
        else
          error('Invalid argument type.');
        end
      end

      % Create partition
      switch partitionType
        case 'KFold'
          obj = obj.createStratifiedKFold(group, K, stratifyOption);
        case 'HoldOut'
          obj = obj.createStratifiedHoldOut(group, P, stratifyOption);
      end
    end

    function obj = createCustomPartition(obj, testSets)
      % testSets is either an Nobsx1 index vector, or an NobsxNset boolean matrix
      assert(isnumeric(testSets)||islogical(testSets));
      if isnumeric(testSets)
        tmp = testSets;
        utmp = unique(tmp);
        testSets = cell(1, numel(utmp));
        for i = 1:numel(testSets)
          testSets{i} = find(tmp==utmp(i));
        end
      else
        tmp = testSets;
        testSets = cell(1, size(testSets,2));
        for i = 1:numel(testSets)
          testSets{i} = find(tmp(:,i));
        end
      end
      
      obj.TestSets        = testSets;
      obj.NumObservations = max(cellfun(@max, testSets));
      obj.NumTestSets     = length(testSets);
      obj.IsCustom        = true;
      obj.Type            = 'CustomPartition';
      obj.TrainSize = obj.NumObservations - cellfun(@length, testSets);
      obj.TestSize = cellfun(@length, testSets);
    end

    function obj = createKFold(obj, N, K)
      % Create K-Fold partition for N observations
      obj.Type = 'KFold';
      obj.NumTestSets = K;
      indices = randperm(N);
      foldSizes = floor(N / K) * ones(1, K);
      foldSizes(1:N - sum(foldSizes)) = foldSizes(1:N - sum(foldSizes)) + 1;
      obj.TestSets = mat2cell(indices, 1, foldSizes);
      obj.TrainSize = N - foldSizes;
      obj.TestSize = foldSizes;
    end

    function obj = createHoldOut(obj, N, P)
      % Create Hold-Out partition
      obj.Type = 'HoldOut';
      obj.NumTestSets = 1;
      if P < 1
        testSize = round(P * N);
      else
        testSize = P;
      end
      indices = randperm(N);
      obj.TestSets = {indices(1:testSize)};
      obj.TrainSize = N - testSize;
      obj.TestSize = testSize;
    end

    function obj = createLeaveOut(obj, N)
      % Create Leave-One-Out partition
      obj.Type = 'LeaveOut';
      obj.NumTestSets = N;
      obj.TestSets = num2cell(1:N);
      obj.TrainSize = N - 1;
      obj.TestSize = 1;
    end

    function obj = createResubstitution(obj, N)
      % Create Resubstitution partition
      obj.Type = 'Resubstitution';
      obj.NumTestSets = 1;
      obj.TestSets = {1:N};
      obj.TrainSize = N;
      obj.TestSize = N;
    end

    function obj = createStratifiedKFold(obj, group, K, stratifyOption)
      % Create K-Fold partition for an array of length(group) observations,
      % with optional stratification

      N = numel(group);
      foldSizes = floor(N / K) * ones(1, K);
      foldSizes(1:N - sum(foldSizes)) = foldSizes(1:N - sum(foldSizes)) + 1;

      % Create Stratified K-Fold partition
      if stratifyOption
        obj.Type     = 'StratifiedKFold';
        obj.TestSets = cell(1, K);

        classes = unique(group);
        numClasses = length(classes);
        classCounts = histcounts(categorical(group), categories(categorical(classes)));

        minNperFold = floor(classCounts/K);
        if any(minNperFold==0)
          warning('not all classes are represented in each of the test folds');
        end

        classIndices = cell(1,K);
        for c = 1:numClasses
          if isnumeric(classes)
            ix = find(group==classes(c));
          else
            ix = find(strcmp(group, classes{c}));
          end
          classIndices{c} = ix(randperm(numel(ix)));
        end

        for c = 1:numClasses
          % Assign to folds what can already be assigned
          for k = 1:K
            endIdx = minNperFold(c);
            obj.TestSets{k} = [obj.TestSets{k}; classIndices{c}(1:endIdx)];
            classIndices{c}(1:endIdx) = [];
          end
        end

        % Now assign the rest
        Nleft = cellfun(@length, classIndices);
        [srt, ix] = sort(Nleft, 'descend');
        classIndices = classIndices(ix);
        for c = 1:numel(classIndices)
          candidates = find(foldSizes - cellfun(@length, obj.TestSets) > 0);
          if numel(candidates) >= srt(c)
            list = randperm(numel(candidates), srt(c));
          else
            list = randi([1 numel(candidates)], [1 srt(c)]);
          end
          for k = list
            obj.TestSets{candidates(k)} = [obj.TestSets{candidates(k)}; classIndices{c}(1)];
            classIndices{c}(1) = [];
          end
        end

        obj.NumTestSets = K;
        obj.TrainSize = obj.NumObservations - cellfun(@length, obj.TestSets);
        obj.TestSize = cellfun(@length, obj.TestSets);
      else
        obj = obj.createKFold(obj.NumObservations, K);
      end
    end

    function obj = createStratifiedHoldOut(obj, group, P, stratifyOption)
      error('this has not yet been tested, and JM does not trust this, given that Mistral messed up the rest as well');

      % Create Stratified Hold-Out partition
      if stratifyOption
        obj.Type = 'StratifiedHoldOut';
        classes = unique(group);
        numClasses = length(classes);
        classCounts = histcounts(categorical(group), categories(categorical(classes)));

        if P < 1
          testSize = round(P * obj.NumObservations);
        else
          testSize = P;
        end

        testIndices = [];
        for c = 1:numClasses
          if isnumeric(classes)
            classIndices = find(group==classes(c));
          else
            classIndices = find(strcmp(group, classes{c}));
          end
          numTestFromClass = round(testSize * (classCounts(c) / obj.NumObservations));
          testIndices = [testIndices; randsample(classIndices, numTestFromClass)];
        end
        obj.TestSets = {testIndices};
        obj.NumTestSets = 1;
        obj.TrainSize = obj.NumObservations - length(testIndices);
        obj.TestSize = length(testIndices);
      else
        obj = obj.createHoldOut(obj.NumObservations, P);
      end
    end

    function testIndices = test(obj, fold)
      % Return test set indices for a given fold
      if fold > obj.NumTestSets
        error('Fold index exceeds the number of test sets.');
      end
      testIndices = false(1, obj.NumObservations);
      testIndices(obj.TestSets{fold}) = true;
    end

    function trainIndices = training(obj, fold)
      % Return training set indices for a given fold
      if fold > obj.NumTestSets
        error('Fold index exceeds the number of test sets.');
      end
      trainIndices = true(1, obj.NumObservations);
      trainIndices(obj.TestSets{fold}) = false;
    end

    function obj = repartition(obj)
      % Rerandomize the partition
      switch obj.Type
        case 'KFold'
          obj = obj.createKFold(obj.NumObservations, obj.NumTestSets);
        case 'HoldOut'
          P = obj.TestSize / obj.NumObservations;
          obj = obj.createHoldOut(obj.NumObservations, P);
        case 'StratifiedKFold'
          obj = obj.createStratifiedKFold(obj.Group, obj.NumTestSets, true);
        case 'StratifiedHoldOut'
          P = obj.TestSize / obj.NumObservations;
          obj = obj.createStratifiedHoldOut(obj.Group, P, true);
        case 'LeaveOut'
          obj = obj.createLeaveOut(obj.NumObservations);
        case 'Resubstitution'
          obj = obj.createResubstitution(obj.NumObservations);
        case 'CustomPartition'
          error('Cannot rerandomize a custom partition.');
      end
    end
  end
end