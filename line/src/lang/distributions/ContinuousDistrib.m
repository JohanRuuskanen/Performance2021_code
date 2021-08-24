classdef ContinuousDistrib < Distrib
    % An abstract class for continuous distributions
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    methods (Hidden)
        function self = ContinuousDistrib(name, numParam, support)
            % SELF = CONTINUOUSDISTRIB(NAME, NUMPARAM, SUPPORT)
            
            % Construct a continuous distribution from name, number of
            % parameters, and range
            self@Distrib(name,numParam,support);
        end
    end
    
end
