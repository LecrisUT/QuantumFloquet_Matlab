classdef rect
    properties
        x
        y
    end
    properties (Dependent)
        w
        h
    end
    methods
        function obj=rect(x,y,Args)
            arguments
                x   (1,2)   double  {mustBeVector}
                y   (1,2)   double  {mustBeVector}
                Args.width  (1,1)   double  {mustBePositive}
                Args.height (1,1)   double  {mustBePositive}
            end
            if (x(2) < x(1) || y(2) < y(1))
                error('Coordinates input are not ordered properly');
            end
            if isfield(Args,'width')
                x(2) = x(1) + Args.width;
            end
            if isfield(Args,'height')
                y(2) = y(1) + Args.height;
            end
            obj.x = x;
            obj.y = y;
        end
        function w=get.w(obj)
            w = obj.x(2) - obj.x(1);
        end
        function h=get.h(obj)
            h = obj.y(2) - obj.y(1);
        end
        function mat=double(obj)
            mat = [obj.x(1) obj.y(1) obj.w obj.h];
        end
        function obj=normalize(obj,inp)
            arguments
                obj
                inp
            end
            % NORMALIZE normalize size to annother object
            % Currently only implemented to normalize to axis position

            switch class(inp)
                case 'matlab.graphics.axis.Axes'
                    % Rename for clarity
                    ax = inp;
                    % Name position variables for clarity
                    pos = ax.Position;
                    px = pos(1);
                    py = pos(2);
                    width = pos(3);
                    height = pos(4);
        
                    switch ax.XScale
                        case 'linear'
                            dx = diff(ax.XLim);
                            x0 = ax.XLim(1);
                            tx = obj.x;
                        case 'log'
                            dx = diff(log10(ax.XLim));
                            x0 = log10(ax.XLim(1));
                            tx = log10(obj.x);
                    end
                    switch ax.YScale
                        case 'linear'
                            dy = diff(ax.YLim);
                            y0 = ax.YLim(1);
                            ty = obj.y;
                        case 'log'
                            dy = diff(log10(ax.YLim));
                            y0 = log10(ax.YLim(1));
                            ty = log10(obj.y);
                    end
                    obj.x = width * ((tx-x0)./dx) + px;
                    obj.y = height * ((ty-y0)./dy) + py;
                otherwise
                    error('Not implemented');
            end
        end
    end
end