classdef GenericCalcIterator < Calc.baseCalcIterator
    % GenericCalcIterator Generic iterator interface
    % Generic implementation of an iterator via function handle. Can not use
    % anonymous functions due to syntax limitations. Use local functions
    % defined in calculation script instead
    %
    % See also baseCalcIterator, function

    properties (SetAccess=protected)
        size
        done
        % data - Data used in the iterator
        % Stored as a vector array
        data        (1,:)
        % updatefcn - Handle of the update function
        % Called on each execution of next via updatefcn(obj)
        % See also function, Calc.GenericCalcIterator.next
        updatefcn
    end
    methods
        function obj = GenericCalcIterator(object,updatefcn,Args)
            arguments
                object
                updatefcn
                Args.data
                Args.size
            end
            % GenericCalcIterator Constructor
            %
            % Syntax:
            %   obj = GenericCalcIterator(object,updatefcn)
            %   [___] = GenericCalcIterator(___,Name,Value)
            % 
            % Description:
            %   obj = GenericCalcIterator(object,updatedfcn) Create an iterator with
            %   updatefcn as the updating function
            %   [___] = GenericCalcIterator(___,Name,Value) specifies options using
            %   name-value arguments in addition to any of the input arguments in
            %   previous syntaxes.
            %
            % Inputs:
            %   object - Initial iterator object
            %   updatefcn - Updating function
            %   Name-Value pairs
            %
            % Outputs:
            %   obj - Iterator object
            %
            % Name-value arguments:
            %   data - Data stored in the iterator
            %   size - [length(data)] Size or the iterator
            %
            % See also GenericCalcIterator, baseCalcIterator,
            % Calc.GenericCalcIterator.next

            obj@Calc.baseCalcIterator(object);
            obj.updatefcn = updatefcn;
            if isfield(Args,'data')
                obj.data = Args.data;
            end
            if isfield(Args,'size')
                obj.size = Args.size;
            else
                obj.size = length(obj.data);
            end
        end
        function val = get.done(obj)
            if obj.ind < obj.size
                val = false;
            elseif obj.ind == obj.size
                val = true;
            else
                error('Out of bounds');
            end
        end
        function next(obj)
            if obj.done
                error('Out of bounds');
            end
            obj.ind = obj.ind + 1;
            obj.updatefcn(obj);
        end
        function reset(obj)
            obj.ind = 0;
        end
    end
end