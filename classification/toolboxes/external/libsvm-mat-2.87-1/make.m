% This make.m is used under Windows

mex -O -c svm.cpp
mex -O -c svm_model_matlab.c
mex -O svmtrain.c svm.obj svm_model_matlab.obj
mex -O svmpredict.c svm.obj svm_model_matlab.obj
mex -O read_sparse.c
