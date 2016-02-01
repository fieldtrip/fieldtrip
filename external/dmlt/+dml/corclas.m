classdef corclas < dml.method
% CORCLAS template matching correlation based classifier.
%
%   DESCRIPTION
%   This classifier makes a template for each class by averaging over
%   instances of that class. A test sample is classifier as the class which
%   corresponds to the template with which the sample has the highest
%   correlation.
%
%   EXAMPLE
%   X = rand(10,20); Y = [1 1 1 1 1 2 2 2 2 2]';
%   m = dml.corclas;
%   m = m.train(X,Y);
%   Z = m.test(X);
%
%   REFERENCE
%   Distributed and overlapping representations of faces and objects in
%   ventral temporal cortex by J. V. Haxby et al.
%
%   DEVELOPER
%   Ali Bahramisharif (ali@cs.ru.nl)

    properties
        template  % the template of each class
    end

    methods

        function obj = corclas(varargin)

            obj = obj@dml.method(varargin{:});

        end

        function obj = train(obj,X,Y)

            nvars = size(X,2);
            nclasses = numel(unique(Y(:,1)));
            obj.template=zeros(nclasses,nvars);
            for i=1:nclasses
                obj.template(i,:)=mean(X(Y==i,:),1);
            end
            
        end

        function Y = test(obj,X)
          
            o=zeros(size(obj.template,1),size(X,1));
            for i=1:size(obj.template,1)
                for j=1:size(X,1)
                    o(i,j)=corr(obj.template(i,:)',X(j,:)');
                end
            end
            Y=o';

        end

        function m = model(obj)
            % MODEL returns the following parameter:
            %
            % m.template

            m.template=obj.template;
            
        end
    end
end
