function tra=fasttransformplanar(data)
% This function is written for fast planar gradient transform using sincos
% method
% The output is the transformation matrix
% Copyright 2009, Ali Bahramisharif

cfgpl.neighbourdist = 4;
montage = megplanar_sincos(cfgpl, data.grad);
chidxO = match_str(montage.labelorg, data.cfg.channel);
chidxN=cat(1,chidxO,chidxO+273);
tra=montage.tra(chidxN,chidxO)';