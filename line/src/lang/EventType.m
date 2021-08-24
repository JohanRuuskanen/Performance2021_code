classdef EventType < Copyable
    % Types of events 
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    % event major classification
    properties (Constant)
        ID_INIT = -1; % model is initialized (t=0)
        ID_LOCAL = 0;
        ID_ARV = 1; % job arrival
        ID_DEP = 2; % job departure
        ID_PHASE = 3; % service advances to next phase, without departure
        ID_READ = 4; % read cache item
        ID_STAGE = 5; % random environment stage change
        
        INIT = categorical("INIT"); 
        LOCAL = categorical("LOCAL");
        ARV = categorical("ARV");
        DEP = categorical("DEP");
        PHASE = categorical("PHASE");
        READ = categorical("READ");
        STAGE = categorical("STAGE");
    end
    
    methods(Static)
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            
            switch type
                case EventType.ARV
                    text = 'ARV';
                case EventType.DEP
                    text = 'DEP';
                case EventType.PHASE
                    text = 'PHASE';
                case EventType.READ
                    text = 'READ';
                case EventType.STAGE
                    text = 'STAGE';
            end
        end
    end
    
end
