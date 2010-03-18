% NMLT is a Matlab machine learning toolbox
% that is specifically tailored to offline and online
% single-trial analysis in cognitive neuroscience.
% 
% This code requires at least Matlab distribution 7.6.0.324 (R2008a)
% 
% Copyrights (C) 2008
% Marcel van Gerven
% Intelligent Systems (http://www.ru.is/ml)
% FC Donders Centre (http://www.ru.nl/fcdonders/)
% Radboud University Nijmegen, The Netherlands
% 
% Most functions in this toolbox are licensed under the GNU General
% Public license (GPL), see http://www.gnu.org for details.
% 
% Unauthorised copying and distribution of functions that are not
% explicitely covered by the GPL is not allowed!
% 
% The functions in this toolbox are copyrighted by their authors:
% 
% Marcel van Gerven
% Pawel Herman
% Ali Bahramisharif
% Jason Farquhar
% Adriana Birlutiu
% Tom Heskes 
%
% Donders Institute for Brain, Cognition and Behaviour
% 
% The toolbox depends on functions from other toolboxes to
% do some of the actual work. These other toolboxes on which
% the framework depends are copyrighted by their respective authors,
% see each individual matlab file for the details.
% 
% gpml-matlab:
% http://www.gaussianprocess.org/gpml/code/matlab/doc/
% 
% L1General:
% http://people.cs.ubc.ca/~schmidtm/Software/L1General/L1General.html
% 
% libsvm:
% libsvm toolbox (http://www.csie.ntu.edu.tw/~cjlin/libsvm/)
% matlab interface (http://www.csie.ntu.edu.tw/~cjlin/libsvm/#matlab)
%
% SLR:
% http://www.cns.atr.jp/~oyamashi/SLR_WEB.html 
%
% hugin (not included in this distribution):
% http://www.hugin.com/Products_Services/Products/Demo/Lite/
%
% For various examples about the functionality of this toolbox, please check
% out the examples folder or use the help function to the methods.
% 
% SEE ALSO:
% base
% validators
% preprocessors
% featureselectors
% classifiers
% regressors
% reconstructors
% 
% TO DO
% 
% - rewrite inner crossvalidation for parameter selection
% - check EM code for mixture models
% - add simple gui
% - integrate utility functions in relevant classes
% - test dynamic_clf and also make dynamic_regr/static_regr
% - repeated calls to train should allow for online learning as much as possible
% - rethink significance testing and also implement for regressor
% - return standard error of mean as second argument of evaluate
% - rewrite feature selection stuff and add searchlight (also check blogreg)
% - better handling of checking of datasets for each method
% - profiling
% - rewrite handling of timeseries...
% - allow arbitrary labels for the classifiers
% - getmodel add descriptions and make sure they can be remapped (e.g.,
%   featureselector)
% - blogreg blinreg return used scale when using multiple updates
% - check imputation in dataset class using e.g. knnimpute
% - clean up metric and significance stuff and split up into clf, regress,
%   and reconstruct
% - what is positive and negative class? be consistent over datasets wrt
%   getmodel
% - test validator.getmodel for transfer learners
% - automate recompilation of mex files
% - create visualizer class
% - IMPLEMENT BLINREG / BLINREG_TRANSFER
%   reshape getmodel whenever possible
% - use separate check_data procedure in mva?
% - check dynamic_classifier and static_classifier code and derivates
% - throw out redundant svm classifiers; check speed and accuracy
% - make one_against_one/rest estimate and map functions suitable for
%   transfer_learning; currently not supported
% - rewrite one_against_one, combiner etc such that transfer learning is
%   used on the separate datasets. This simplifies code and allows for more complex designs 
% - refactor using: perl -p -i -e 's/clfmethod/mvmethod/g' `grep -lr -e 'clfmethod' *`
% - how to handle multiple regression vs regressors/classifiers that return
%   mean + error estimates? see e.g. the kalman filter
% - clean up nclasses field
% - add wilcoxon test; add multiple outputs for evaluate/significance
% - feature selection methods should return sparse matrices instead of
%   reduced matrices in order to keep dimensions consistent. This also allows
%   models for different folds to be combined (currently fails due to
%   incompatible dimensions). This requires that all methods can deal with
%   sparse inputs where the zero inputs are disregarded
% - crossvalidator 276: roundoff error changes the actual number of used
%   samples; can go wrong with skewed unique samples; also for ten-fold