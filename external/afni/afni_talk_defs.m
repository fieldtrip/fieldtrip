function [comm] = afni_talk_defs(comm)
  
   comm.con = -1;
   comm.sockcon = -1;
   comm.socknum = 53211;
   comm.dbg = 0;
   comm.mm = {};
   comm.tt = {};
   
   return;
