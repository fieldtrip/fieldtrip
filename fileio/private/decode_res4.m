function R4F = decode_res4(chunk)

% decode the res4 file content
if 0
  % using  READ_CTF_RES4 -- this does not produce a proper .grad structure
  % TODO: remove this code, and possibly read_ctf_res4 as well
  tmp_name = tempname;
  F = fopen(tmp_name, 'wb');
  fwrite(F, chunk, 'uint8');
  fclose(F);
  R4F = read_ctf_res4(tmp_name);
  delete(tmp_name);
else
  % using FT_READ_HEADER recursively, and then readCTFds in the second call
  % this will also call ctf2grad and yield correct gradiometer information
  tmp_name = tempname;
  [dirname, fname] = fileparts(tmp_name);
  res4fn = [tmp_name '.ds/' fname '.res4'];
  meg4fn = [tmp_name '.ds/' fname '.meg4'];
  dsname = [tmp_name '.ds'];
  
  mkdir(dsname);
  
  F = fopen(res4fn, 'wb');
  fwrite(F, chunk, 'uint8');
  fclose(F);
  
  F = fopen(meg4fn, 'wb');
  fwrite(F, 'MEG42CP');
  fclose(F);
  
  R4F = ft_read_header(dsname, 'coordsys', 'dewar');
  
  delete(res4fn);
  delete(meg4fn);
  rmdir(dsname);
end

