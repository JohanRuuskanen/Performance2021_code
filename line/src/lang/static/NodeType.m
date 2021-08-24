classdef (Sealed) NodeType
    % Enumeration of node types.
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        Transition = 9;
        Place = 8;
        Fork = 7;
        Router = 6;
        Cache = 5;
        Logger = 4;
        ClassSwitch = 3;
        Delay = 2;
        Source = 1;
        Queue = 0;
        Sink = -1;
        Join = -2;
    end
    
    methods (Static)
        function bool = isStation(nodetype)
            % BOOL = ISSTATION(NODETYPE)
            
            bool = (nodetype == NodeType.Source | nodetype == NodeType.Delay | nodetype == NodeType.Queue | nodetype == NodeType.Join | nodetype == NodeType.Place);
        end
        function bool = isStateful(nodetype)
            % BOOL = ISSTATEFUL(NODETYPE)            
            bool = (nodetype == NodeType.Source | nodetype == NodeType.Delay | nodetype == NodeType.Queue | nodetype == NodeType.Cache | nodetype == NodeType.Join | nodetype == NodeType.Router | nodetype == NodeType.Place);
        end
    end
    
end
