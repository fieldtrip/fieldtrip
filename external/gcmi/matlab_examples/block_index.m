function idx = block_index(bi, blocklen, siglen)

start = (bi-1)*blocklen + 1;
finish = min(start+blocklen-1, siglen);

idx = start:finish;