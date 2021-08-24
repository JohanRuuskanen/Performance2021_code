classdef (Sealed) TimingStrategy
    % Enumeration of timing polcies in petri nets transitions.
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        Timed = 0;
        Immediate = 1;
    end
    
    methods (Static)
        
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            
            switch type
                case TimingStrategy.Timed
                    text = 'Timed Transaction';
                case TimingStrategy.Immediate
                    text = 'Immediate Transaction';
            end
        end
    end
end

