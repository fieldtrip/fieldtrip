classdef KTFigure
    properties (Hidden, SetAccess = private, GetAccess = private)
        figureHandle
    end
    methods (Hidden)
        function obj = KTFigure
            obj.figureHandle = figure;
            set(obj.figureHandle, 'Visible', 'off');
        end
        function screenSize = getScreenSize(~)
            screenSize = get(0, 'ScreenSize');
        end
    end
    methods
        function handle = getHandle(obj)
            handle = obj.figureHandle;
        end
        function EnlargeFigure(obj)
            screenSize = obj.getScreenSize;
            set(gcf, 'position', [screenSize(3:4).*[1.7, 1]*0.1 screenSize(3:4).*[1, 1.3]*0.6]);
        end
        function EnlargeFigureWide(obj)
            screenSize = obj.getScreenSize;
            set(gcf, 'position', [screenSize(3:4).*[0.3, 1]*0.1 screenSize(3:4).*[1.55, 1.3]*0.6]);
        end
        function HideFigure(obj)
            set(obj.figureHandle, 'Visible', 'off');
        end
        function ShowFigure(obj)
            set(obj.figureHandle, 'Visible', 'on');
        end
        function MakeBackgroundWhite(obj)
            set(obj.figureHandle, 'color', [1 1 1]);
        end
        function MakeBackgroundGray(obj)
            set(obj.figureHandle, 'color', [0.5 0.5 0.5]);
        end
        function MakeBackgroundBlack(obj)
            set(obj.figureHandle, 'color', [0 0 0]);
        end
        function SetActive(obj)
            figure(obj.figureHandle);
        end
        function CloseFigure(obj)
            close(obj.figureHandle);
        end
    end
end
