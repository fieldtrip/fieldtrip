function newmodel = submodel(model,dims,mask,req)
% SUBMODEL retrieves a submodel from the full model:
% 
% newmodel is the new model in the requested form
%
% model = the data as a vector
% dims  = original dimensions of the model
% mask  = linear mask indicating which elements are specified
% req   = vector specifying which elements are required per dimension
%         nan indicates the full range for that dimension

  newmodel = [];

  didx = isnan(req);
  
  newdims = dims(didx);

  if isempty(newdims)
    error('empty dimensions');
  end

  if length(newdims) == 1
    newdims = [1 newdims];
  end
 
  orig = ind2subv(dims,mask);

  if all(didx)
    sameidx = 1:size(orig,1);
  else
    sameidx = find(ismember(orig(:,~didx),req(~didx),'rows'));
  end

  if ~isempty(sameidx)

    % non-sparse new model
    newmodel = zeros(newdims);
    
    for j=1:length(sameidx)
      newmodel(subv2ind(newdims,orig(sameidx(j),didx))) = model(sameidx(j));
    end
  end
  
end





