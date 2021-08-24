classdef (Sealed) DropStrategy
    % Enumeration of drop policies in stations.
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        InfiniteBuffer = -1;
        Drop = 1;
        BlockingAfterService = 2;
    end
    
    methods (Static)
        
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            
            switch type
                case DropStrategy.InfiniteBuffer
                    text = 'waiting queue';
                case DropStrategy.Drop
                    text = 'drop';
                case DropStrategy.BlockingAfterService
                    text = 'BAS blocking';
            end
        end
    end

end

