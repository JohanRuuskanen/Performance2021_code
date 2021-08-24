classdef Processor  < Host
    % A hardware server in a LayeredNetwork.
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    
    methods
        %public methods, including constructor
        
        %constructor
        function self = Processor(model, name, multiplicity, scheduling, quantum, speedFactor)
            % OBJ = PROCESSOR(MODEL, NAME, MULTIPLICITY, SCHEDULING, QUANTUM, SPEEDFACTOR)
            if ~exist('name','var')
                line_error(mfilename,'Constructor requires to specify at least a name.');
            end
            
            if ~exist('multiplicity','var')
                multiplicity = 1;
            end
            if ~exist('scheduling','var')
                scheduling = SchedStrategy.PS;
            end
            if ~exist('quantum','var')
                quantum = 0.001;
            end
            if ~exist('speedFactor','var')
                speedFactor = 1;
            end            
            self@Host(model, name, multiplicity, scheduling, quantum, speedFactor)
        end
        
    end
    
end
