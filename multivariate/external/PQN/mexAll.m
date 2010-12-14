fprintf('Compiling minFunc files...\n');
mex minFunc/lbfgsC.c
fprintf('Compiling UGM files...\n');
mex -IUGM/mex UGM/mex/UGM_makeNodePotentialsC.c
mex -IUGM/mex UGM/mex/UGM_makeEdgePotentialsC.c
mex -IUGM/mex UGM/mex/UGM_PseudoLossC.c
mex -IUGM/mex UGM/mex/UGM_Infer_ExactC.c
mex -IUGM/mex UGM/mex/UGM_Loss_subC.c
mex -IUGM/mex UGM/mex/UGM_updateGradientC.c
mex -IUGM/mex UGM/mex/UGM_Sample_GibbsC.c
mex -IUGM/mex UGM/mex/UGM_MRFLoss_subC.c
fprintf('Compiling projection files...\n');
mex project/projectRandom2C.c
mex -Iproject project/projectBlockL1.c project/oneProjectorCore.c project/heap.c
mex -Iproject project/projectBlockL2.c
fprintf('Compiling CRF files (requires a good compilier)...\n\n');
mex crfChain/crfChain_lossL2C.c