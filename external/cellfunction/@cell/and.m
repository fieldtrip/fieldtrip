function out = and(in1, in2)

out = cellfun(@and,in1,in2,'uniformoutput',false);
