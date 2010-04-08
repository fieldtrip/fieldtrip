
function svm_out = svm_eval(data,sv,weights,bias,kernel,kerparam)

 K  = calc_kernel(data,sv,kernel,kerparam);
 svm_out = K * weights + bias;