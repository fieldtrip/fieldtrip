function gui_help(which_help)
%
% Used by FASTICAG

% All the help texts and title used by GUI are stored here. 
% Make changes here.
% Also displays the helpwindow with the selected text

% @(#)$Id: gui_help.m,v 1.6 2005/10/19 13:05:34 jarmo Exp $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch which_help
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'pcamat'

  helptitle = 'FastICA: Reduce dimension';
  helptext=[ ...
    'You may reduce the dimension of the data by selecting only the     '
    'subspace corresponding to certain eigenvalues of the covariance    '
    'matrix of the data. Give the indices of the first and last         '
    'eigenvalues (sorted in descending order) to be included (all       '
    'eigenvalues in between will be included as well).  The eigenvalues '
    'and their indices can be seen in the graphical plot now on the     '
    'screen. The heights of the bars give the eigenvalues, with indices '
    'below.                                                             '
    '                                                                   '
    'For example, give ''1'' and ''n'' if you want to reduce the dimension  '
    'to n by principal component analysis, which means discarding the   '
    'subspaces corresponding to the smallest eigenvalues. Such a        '
    'dimension reduction may reduce noise and improve the performance of'
    'ICA.                                                               '];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'gui_cb_about'

  helptitle='About FastICA';
  helptext =[ ...
    'FastICA for Matlab 7.x and 6.x                                            '
    'Version 2.5, October 19 2005                                              '
    'Copyright (c) Hugo Gävert, Jarmo Hurri, Jaakko Särelä, and Aapo Hyvärinen.'
    '                                                                          '
    'For more information please see:                                          '
    'http://www.cis.hut.fi/projects/ica/fastica/                               '];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'gui_cb_help'

  helptitle='FastICA GUI';
  helptext = [...
    'Basic function:                                                            '
    '                                                                           '
    '- Click LOAD DATA and give the name of the variable that contains          '
    '  the data.                                                                '
    '                                                                           '
    '- Click DO ICA to perform the analysis.                                    '
    '                                                                           '
    '- Click SAVE RESULTS to store the results for future use.                  '
    '                                                                           '
    'Options:                                                                   '
    '                                                                           '
    'If the input matrix contains the signals as column vectors instead of      '
    'row vectors, click on TRANSPOSE to transpose the data matrix.              '
    '                                                                           '
    'Click on PLOT DATA to see the data as 1-D time signals.                    '
    '                                                                           '
    'Clicking REDUCE DIM gives you a graphical plot of the eigenvalue           '
    'structure of the covariance matrix of the data. You can then reduce        '
    'the dimension of the data by retaining only the subspaces corresponding to '
    'the largest (or smallest) eigenvalues (i.e. variances). To undo this       '
    'operation click ORIGINAL DIM. You can plot the whitened (preprocessed      '
    'data) by PLOT WHITENED.                                                    '
    '                                                                           '
    'Click on DO ICA to perform independent component analysis.                 '
    'Clicking on PLOT ICS has the same effect, except that DO ICA forces        '
    'recomputation of ICA.                                                      '
    '                                                                           '
    'You can choose the decorrelation approach by the ''Approach'' drop-down menu:'
    'deflation means that the independent components are estimated              '
    'one-by-one, whereas in the symmetric approach they are estimated in        '
    'parallel. You can now choose the number of independent components to be    '
    'estimated in both deflation and symmetric approaches.                      '
    '                                                                           '
    'You have a choice of three nonlinearities:                                 '
    '                                                                           '
    '''pow3'' (default) :  g(u)=u^3                                               '
    '''tanh''           :  g(u)=tanh(u)                                           '
    '''gauss''          :  g(u)=u*exp(-u^2/2)                                     '
    '''skew''           :  g(u)=u^2                                               '
    '                                                                           '
    'For example, you could choose approach=''symmetric'' and nonlinearity=''tanh'' '
    'to perform maximum likelihood ICA estimation for supergaussian data.       '
    '                                                                           '
    'If the algorithm does not seem to converge, you can use the stabilized     '
    'version of the fixed-point algorithm. To use the stabilized version,       '
    'choose ''on'' from the drop-down menu ''Stabilization''.                       '
    'If you have specified a value less than 1 for the parameter ''mu'' from      '
    'the ''Advanced Options'' menu then the ''Stabilization'' drop-down menu is     '
    'not active. This is because if the parameter ''mu'' is less than 1 then the  '
    'program will use the stabilized code. Please see the help for              '
    'Advanced Options for more information about stabilization.                 '
    '                                                                           '
    'The ADVANCED OPTIONS menu has its own HELP button.                         '
    '                                                                           '
    'During computations, an INTERRUPT button appears. Clicking the button      '
    'interrupts the computations.                                               '];
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'gui_advc'

  helptitle='FastICA GUI: Advanced options';
  helptext = [...
    'Advanced options:                                                     '
    '                                                                      '
    'In some cases, it may be desired to improve the statistical           '
    'performance of the algorithm by using a fine-tuning procedure. This   '
    'means that after (initial) convergence, the algorithm is run          '
    '(possibly) with a different nonlinearity, and using a smaller step    '
    'size and the stabilized version of the fixed-point algorithm.         '
    'You can specify the nonlinearity that will then be used with the      '
    'parameter ''Finetune''. If you set the the finetuning to ''off'', then the'
    'fine-tuning won''t be done.                                            '
    '                                                                      '
    'You can also fine-tune the nonlinearities used in the fixed-point     '
    'algorithm.                                                            ' 
    'The nonlinearities tanh and gauss contain parameters a1 and a2, so    '
    'that the nonlinearities are in fact defined as:                       '
    '''tanh''           :  g(u)=tanh(a1*u)                                   '
    '''gauss''          :  g(u)=u*exp(-a2*u^2/2)                             '
    'The default values of a1 and a2 are 1, in which case they effectively '
    'disappear from the definitions.                                       '
    '                                                                      '
    'If the algorithm does not seem to converge, you can use the stabilized'
    'version of the fixed-point algorithm. There are two ways of doing     '
    'this. The first one is to explicitly specify the value of the step    '
    'size parameter ''mu''. The default value is 1. Choosing a value that    '
    'is smaller than 1 implies that the computations are made using the    '
    'stabilized fixed-point algorithm. The second way to use the stabilized'
    'version is simpler: choose ''on'' in the drop-down menu ''stabilization'' '
    '(on the main menu page). Then the value of mu will be changed         '
    'automatically during the ICA calculations. If the program senses that '
    'the algorithm is stuck between two points, it will halve the value of '
    'mu (.5 * mu) for duration of one round. (This is called a ''stroke.''). '
    'Also if there is no convergence before half of the maximum number of  '
    'iterations has been reached then the mu will be halved for the rest   ' 
    'of the rounds.                                                        '
    '                                                                      '
    'The parameter ''epsilon'' is used to decide if the algorithm has        '
    'converged. A larger epsilon makes the convergence test less strict.   '
    'Note that if you use finetuning or stabilization, epsilon may need to '
    'be reduced accordingly.                                               '
    '                                                                      '
    '''Maximum number of iterations'' gives the absolute maximum of          '
    'iterations used in the estimation procedure. In the deflation         '
    'approach, this is iterations per component.                           '
    '                                                                      '
    'You can input the ''Initial state'' of the algorithm, i.e. the initial  '
    'value for A. Choose ''guess'' in the drop-down menu ''Initial state'',    '
    'click on ''Load Initial guess'', and give the name of the variable in   '
    'Matlab workspace that contains the initial value.                     '
    '                                                                      '
    'In the drop-down menu ''Display mode'' you can choose if the results are'
    'plotted during computations.  You may wish to switch this off         '
    'especially if you have lots of data which takes a long time to plot.  '
    '''Iteration between displays'' tells how often the running estimates of '
    'the independent components are plotted: A value of 1 means after every'
    'iteration.                                                            '
    '                                                                      '
    'If the data vector is very long (more than 10 000 points), it may be  '
    'advisable to use only a part of the data at every iteration. The      '
    'option ''Sample size'' allows you to give the proportion (0-1) of the   '
    'data that is used at every step. The sample is chosen randomly at     '
    'every step.                                                           '
    '                                                                      '
    'Click on DEFAULT to return to default values for all advanced options.'
    'You can make the new values take effect without closing the window by '
    'clicking APPLY.                                                       '];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'gui_lc_data'

  helptitle='FastICA GUI: Load data';
  helptext = [...
    'Input the name of the variable in Matlab workspace that contains the  '
    'data. The data must be in a single matrix, each row (or column) giving'
    'the values of one signal. If the signals are in column vectors, click '
    'TRANSPOSE after loading the data to transpose the data matrix.        '
    '                                                                      '
    'If the data is in a file, load it to Matlab workspace first.          '];

case 'gui_lc_guess'

  helptitle='FastICA GUI: Load guess';
  helptext = [...
    'Input the name of the variable in Matlab workspace that contains the'
    'initial value for the mixing matrix A, and click OK. If the initial '
    'value is in a file, load it to Matlab workspace first.              '];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'gui_sc'

  helptitle='FastICA GUI: Save results';
  helptext = [...
    'The results will be saved as variables in Matlab workspace.            '
    'You give a suffix that identifies the variables. For example, if you   '
    'give ''_FASTICA'', the results will be stored in the following variables:'
    '                                                                       '
    'W_FASTICA   : estimate of the separating matrix                        '
    'A_FASTICA   : estimate of the mixing matrix                            '
    'IC_FASTICA  : estimated independent components (row vectors)           '
    '                                                                       '
    'Additional results related to preprocessing:                           '
    'D_FASTICA and E_FASTICA    : give the eigenvalue decomposition of the  '
    '                             covariance matrix                         '
    'whiteningMatrix_FASTICA    : matrix performing whitening and dimension '
    '                             reduction                                 '
    'dewhiteningMatrix_FASTICA  : the pseudoinverse of the whitening matrix '
    'whitesig_FASTICA           : whitened (i.e. preprocessed) signals.     '];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

helpwin(helptext, helptitle);