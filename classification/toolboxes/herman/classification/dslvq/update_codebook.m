function codebook = update_codebook(cods,data,coeff)

data_vec = ones(size(cods,1),1) * data;
codebook = cods + coeff.*(data_vec - cods);