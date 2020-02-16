function v = randCauchy(sz,scale,location)
% Return random numbers from the Cauchy distribution 
% with scale and location
% INPUT
% sz - size of the random numbers to return [1 1]
% scale - scale of the Cauchy  [1]
% location - Location of the Cauchy [0]

if nargin<3
    location =0;
    if nargin<2
        scale =1;
        if nargin <1 
            sz = [1 1];
        end
    end
end

p = rand(sz);
v =location+ scale*tan(pi*(p-0.5)); % Use the inverse CDF to generate random numbers.

