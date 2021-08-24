classdef (Sealed) RoutingStrategy
    % Enumeration of routing strategies
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        ID_RAND = 0;
        ID_PROB = 1;
        ID_RRB = 2;
        ID_JSQ = 3;
        ID_FIRING = 4;
        ID_DISABLED = -1;
        RAND = 'Random';
        RRB = 'RoundRobin';
        PROB = 'Probabilities';
        JSQ = 'JoinShortestQueue';
        DISABLED = 'Disabled';
        FIRING = 'Firing';
    end
    
    methods (Static, Access = public)
        function type = toType(text)
            % TYPE = TOTYPE(TEXT)
            
            switch text
                case 'Random'
                    type = RoutingStrategy.RAND;
                case 'RoundRobin'
                    type = RoutingStrategy.RRB;
                case RoutingStrategy.PROB
                    type = RoutingStrategy.PROB;
                case 'JoinShortestQueue'
                    type = RoutingStrategy.JSQ;
                case 'Disabled'
                    type = RoutingStrategy.DISABLED;
                case 'Firing'
                    type = RoutingStrategy.FIRING;
            end
        end
        
        function feature = toFeature(type)
            % FEATURE = TOFEATURE(TYPE)            
            switch type
                case RoutingStrategy.RAND
                    feature = 'RoutingStrategy_RAND';
                case RoutingStrategy.RRB
                    feature = 'RoutingStrategy_RRB';
                case RoutingStrategy.PROB
                    feature = 'RoutingStrategy_PROB';
                case RoutingStrategy.JSQ
                    feature = 'RoutingStrategy_JSQ';
                case RoutingStrategy.DISABLED
                    feature = 'RoutingStrategy_DISABLED';
                case RoutingStrategy.FIRING
                    feature = 'RoutingStrategy_FIRING';
                case 0 % if unassigned, set it by default to Disabled
                    feature = 'RoutingStrategy_DISABLED';
            end
        end
        
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            
            switch type
                case RoutingStrategy.RAND
                    text = 'Random';
                case RoutingStrategy.RRB
                    text = 'RoundRobin';
                case RoutingStrategy.PROB
                    text = RoutingStrategy.PROB;
                case RoutingStrategy.JSQ
                    text = 'JoinShortestQueue';
                case RoutingStrategy.DISABLED
                    text = 'Disabled';
                case RoutingStrategy.FIRING
                    text = 'Firing';
                case 0 % if unassigned, set it by default to Disabled
                    text = 'Disabled';
            end
        end
    end
    
end

