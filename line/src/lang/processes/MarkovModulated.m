classdef MarkovModulated < PointProcess
    % An abstract class for Markov-modulated processes
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    methods (Hidden)
        %Constructor
        function self = MarkovModulated(name, numParam)
            % SELF = MARKOVMODULATED(NAME, NUMPARAM)
            
            self@PointProcess(name, numParam);
        end
    end
    
    methods
        function X = sample(self, n)
            % X = SAMPLE(N)
            
            if ~exist('n','var'), n = 1; end
            MAP = self.getRepresentation;
            if map_isfeasible(MAP)
                X = map_sample(MAP,n);
            else
                line_error(mfilename,'This process is infeasible (negative rates).');
            end
        end
    end
    
    methods %(Abstract) % implemented with errors for Octave compatibility
        function phases = getNumberOfPhases(self)
            % PHASES = GETNUMBEROFPHASES()
            
            line_error(mfilename,'Line:AbstractMethodCall','An abstract method was called. The function needs to be overridden by a subclass.');
        end
        function MAP = getRepresentation(self)
            % MAP = GETREPRESENTATION()
            
            line_error(mfilename,'Line:AbstractMethodCall','An abstract method was called. The function needs to be overridden by a subclass.');
        end
    end
    
    methods (Static)
        function cx = fit(MEAN, SCV)
            % CX = FIT(MEAN, SCV)
            
            line_error(mfilename,'Line:AbstractMethodCall','An abstract method was called. The function needs to be overridden by a subclass.');
        end
    end
    
end

