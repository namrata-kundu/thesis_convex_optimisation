classdef SingleVariablePLQ
    %SINGLEVARIABLEPLQ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        f
        pieces
        rho
        new_pieces
    end
    
    methods
        function obj = SingleVariablePLQ(func, pieces)
            %SINGLEVARIABLEPLQ Construct an instance of this class
            %   Detailed explanation goes here
            obj.f = func;
            obj.pieces = pieces;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

