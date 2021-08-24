classdef LayeredNetworkElement < Element
    % A generic element of a LayeredNetwork model.
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        HOST = 0;
        PROCESSOR = 0;
        TASK = 1;
        ENTRY = 2;
        ACTIVITY =3;
        CALL = 4;
    end
    
    
    methods
        %Constructor
        function self = LayeredNetworkElement(name)
            % SELF = LAYEREDNETWORKELEMENT(NAME)
            
            self@Element(name);
        end
        
    end
    
end
