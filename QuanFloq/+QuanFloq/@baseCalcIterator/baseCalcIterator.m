classdef baseCalcIterator < handle & matlab.mixin.Heterogeneous
    % baseCalcIterator Base iterator interface
    % Simple iterator interface for adiabatic calculations or other such
    % applications
    %
    % See also baseCalc, GenericCalcIterator
    
    properties (SetAccess=protected)
        % ind - Current index of the iteration
        ind     (1,1)   double  {mustBeInteger} =0
        % object - Current Calc object generated
        object  QuanFloq.baseCalc   {mustBeScalarOrEmpty}
    end
    properties (Abstract, SetAccess=protected)
        % size - Maximum size of the iterator
        size    (1,1)   double  {mustBeInteger}
        % done - Whether or not the iterator has finished
        done    (1,1)   logical
    end
    methods (Abstract)
        % next - Get next object in iterator
        next(obj)
        % reset - Reset the iterator to initial state
        reset(obj)
    end
    methods
        function obj = baseCalcIterator(object)
            arguments
                object
            end
            % baseCalcIterator Constructor
            %
            % Syntax:
            %   obj = baseCalcIterator(object)
            % 
            % Description:
            %   obj = baseCalcIterator(object) Create the base iterator with object as
            %   initial baseCalc object
            %
            % Inputs:
            %   object - Initial iterator object
            %
            % Outputs:
            %   obj - Iterator object
            % 
            % See also baseCalcIterator
            
            obj.object = object;
        end
    end
end