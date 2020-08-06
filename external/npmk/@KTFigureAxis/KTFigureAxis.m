classdef KTFigureAxis
    properties (Access = public)
        linSpaceX = 20;
        linSpaceY = 20;
        lowerLimitX
        lowerLimitY
        upperLimitX
        upperLimitY
        axisHandle
    end
    methods (Hidden)
        function obj = KTFigureAxis(lowerLimitX, upperLimitX, lowerLimitY, upperLimitY, axisHandle)
            switch nargin
                case 0
                    obj.axisHandle = gca;
                    disp('Set the limits.');
                case 4
                    obj.lowerLimitX = lowerLimitX;
                    obj.upperLimitX = upperLimitX;
                    obj.lowerLimitY = lowerLimitY;
                    obj.upperLimitY = upperLimitY;
                    obj.axisHandle  = gca;
                case 5
                    obj.lowerLimitX = lowerLimitX;
                    obj.upperLimitX = upperLimitX;
                    obj.lowerLimitY = lowerLimitY;
                    obj.upperLimitY = upperLimitY;
                    obj.axisHandle  = axisHandle;                  
                otherwise
                    disp('Number of input arguments is invalid.');
            end
        end
    end
    methods
        function obj = set.lowerLimitX(obj, limit)
            obj.lowerLimitX = limit;
        end
        function obj = set.upperLimitX(obj, limit)
            obj.upperLimitX = limit;
        end
        function obj = set.lowerLimitY(obj, limit)
            obj.lowerLimitY = limit;
        end
        function obj = set.upperLimitY(obj, limit)
            obj.upperLimitY = limit;
        end
        function obj = set.axisHandle(obj, limit)
            obj.axisHandle = limit;
        end
        function obj = set.linSpaceX(obj, spacing)
            obj.linSpaceX = spacing;
        end
        function obj = set.linSpaceY(obj, spacing)
            obj.linSpaceY = spacing;
        end
        function setXLabels(obj)
            xLim = get(gca, 'XLim');
            XTickMarkLocations = linspace(xLim(1), xLim(2), obj.linSpaceX);
            XTickMarkLabelsLocations = linspace(obj.lowerLimitX, obj.upperLimitX, obj.linSpaceX);
            for labelIDX = 1:obj.linSpaceX
                XTickMarkLabels{labelIDX} = num2str(roundn(XTickMarkLabelsLocations(labelIDX),1));
            end
            set(obj.axisHandle, ...
                'XTick', XTickMarkLocations,...
                'XTickLabel', XTickMarkLabels);
        end
        function setYLabels(obj)
            yLim = get(gca, 'YLim');
            YTickMarkLocations = linspace(yLim(1), yLim(2), obj.linSpaceY);
            YTickMarkLabelsLocations = linspace(obj.lowerLimitY, obj.upperLimitY, obj.linSpaceY);
            for labelIDX = 1:obj.linSpaceY
                YTickMarkLabels{labelIDX} = num2str(roundn(YTickMarkLabelsLocations(labelIDX),1));
            end
            set(obj.axisHandle, ...
                'YTick', YTickMarkLocations,...
                'YTickLabel', YTickMarkLabels);
        end
    end
end