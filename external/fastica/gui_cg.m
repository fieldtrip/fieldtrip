function gui_cg()
%
% This function is needed by FASTICAG

% This file just removes the global variables
% that are used in FASTICAG from the memory

% @(#)$Id: gui_cg.m,v 1.2 2003/04/05 14:23:57 jarmo Exp $

clear global c_FastICA_appr_strD;
clear global c_FastICA_appr_strV;
clear global c_FastICA_dMod_strD;
clear global c_FastICA_dMod_strV;
clear global c_FastICA_finetune_strD;
clear global c_FastICA_finetune_strV;
clear global c_FastICA_g1_strD;
clear global c_FastICA_g1_strV;
clear global c_FastICA_g2_strD;
clear global c_FastICA_g2_strV;
clear global c_FastICA_iSta_strD;
clear global c_FastICA_iSta_strV;
clear global c_FastICA_stabili_strD;
clear global c_FastICA_stabili_strV;
clear global c_FastICA_verb_strD;
clear global c_FastICA_verb_strV;
clear global g_FastICA_a1;
clear global g_FastICA_a2;
clear global g_FastICA_approach;
clear global g_FastICA_displayIn;
clear global g_FastICA_displayMo;
clear global g_FastICA_epsilon;
clear global g_FastICA_finetune;
clear global g_FastICA_g;
clear global g_FastICA_ica_A;
clear global g_FastICA_ica_W;
clear global g_FastICA_ica_sig;
clear global g_FastICA_initGuess;
clear global g_FastICA_initState;
clear global g_FastICA_interrupt;
clear global g_FastICA_loadType;
clear global g_FastICA_maxFinetune;
clear global g_FastICA_maxNumIte;
clear global g_FastICA_mixedmean;
clear global g_FastICA_mixedsig;
clear global g_FastICA_myy;
clear global g_FastICA_numOfIC;
clear global g_FastICA_pca_D;
clear global g_FastICA_pca_E;
clear global g_FastICA_sampleSize;
clear global g_FastICA_stabilization;
clear global g_FastICA_verbose;
clear global g_FastICA_white_dwm;
clear global g_FastICA_white_sig;
clear global g_FastICA_white_wm;
clear global hb_FastICA_initGuess;
clear global he_FastICA_a1;
clear global he_FastICA_a2;
clear global he_FastICA_displayInterval;
clear global he_FastICA_epsilon;
clear global he_FastICA_file;
clear global he_FastICA_suffix;
clear global he_FastICA_maxFinetune;
clear global he_FastICA_maxIterations;
clear global he_FastICA_myy;
clear global he_FastICA_numOfIC;
clear global he_FastICA_sampleSize;
clear global hf_FastICA_MAIN;
clear global hf_FastICA_AdvOpt;
clear global hf_FastICA_Load;
clear global hf_FastICA_Save;
clear global hpm_FastICA_approach;
clear global hpm_FastICA_displayMode;
clear global hpm_FastICA_finetune;
clear global hpm_FastICA_g;
clear global hpm_FastICA_initState;
clear global hpm_FastICA_stabilization;
clear global hpm_FastICA_verbose;
clear global ht_FastICA_dim;
clear global ht_FastICA_icaStatus;
clear global ht_FastICA_initGuess;
clear global ht_FastICA_mixedStatus;
clear global ht_FastICA_newDim;
clear global ht_FastICA_numOfSamp;
clear global ht_FastICA_whiteStatus;
