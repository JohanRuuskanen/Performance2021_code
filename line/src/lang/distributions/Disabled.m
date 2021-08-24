classdef Disabled < ContinuousDistrib & DiscreteDistrib
    % A distribution that is not configured
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    methods
        function self = Disabled()
            % SELF = DISABLED()
            
            % Constructs a disabled distribution
            self@ContinuousDistrib('Disabled',1,[NaN,NaN]);
            self@DiscreteDistrib('Disabled',1,[NaN,NaN]);
            setParam(self, 1, 'value', NaN, 'java.lang.Double');
        end
    end
    
    methods
        function bool = isContinuous(self)
            % BOOL = ISCONTINUOUS()
            
            % Returns true is the distribution is continuous
            bool = true;
        end
        
        function bool = isDiscrete(self)
            % BOOL = ISDISCRETE()
            
            % Returns true is the distribution is discrete
            bool = true;
        end
        
        function X = sample(self, n)
            % X = SAMPLE(N)
            
            % Get n samples from the distribution
            if ~exist('n','var'), n = 1; end
            X = nan(n,1);
        end
        
        function ex = getMean(self)
            % EX = GETMEAN()
            
            % Get distribution mean
            ex = NaN;
        end
        
        function SCV = getSCV(self)
            % SCV = GETSCV()
            
            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            SCV = NaN;
        end
        
        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)
            
            % Evaluate the cumulative distribution function at t
            % AT T
            
            Ft = NaN;
        end
        
        function p = evalPMF(self, k)
            % P = EVALPMF(K)
            
            % Evaluate the probability mass function at k
            % AT K
            
            p = 0*k;
        end
    end
    
end

