classdef DiscreteDistrib < Distrib
    % An abstract class for continuous distributions
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    methods (Hidden)
        function self = DiscreteDistrib(name, numParam, support)
            % SELF = DISCRETEDISTRIB(NAME, NUMPARAM, SUPPORT)
            
            % Construct a continuous distribution from name, number of
            % parameters, and range
            self@Distrib(name,numParam,support);
        end
    end
    
end
