classdef Immediate < Distrib
    % A distribution with probability mass entirely at zero
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    methods (Hidden)
        %Constructor
        function self = Immediate()
            % SELF = IMMEDIATE()
            
            self@Distrib('Immediate', 0,[0,0]);
        end
    end
    
    methods
        function bool = isDisabled(self)
            % BOOL = ISDISABLED()
            
            bool = false;
        end
        
        function X = sample(self, n)
            % X = SAMPLE(N)
            
            if ~exist('n','var'), n = 1; end
            X = zeros(n,1);
        end
        
        function ex = getMean(self)
            % EX = GETMEAN()
            
            % Get distribution mean
            ex = 0;
        end
        
        function SCV = getSCV(self)
            % SCV = GETSCV()
            
            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            
            
            SCV = 0;
        end
        
        function mu = getMu(self)
            % MU = GETMU()
            
            % Return total outgoing rate from each state
            mu = Distrib.InfRate;
        end
        
        function phi = getPhi(self)
            % PHI = GETPHI()
            
            % Return the probability that a transition out of a state is
            % absorbing
            phi = 1.0;
        end
        
        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)
            
            % Evaluate the cumulative distribution function at t
            % AT T
            
            Ft = 1;
        end
        
        function L = evalLST(self, s)
            % L = EVALST(S)
            
            % Evaluate the Laplace-Stieltjes transform of the distribution function at t            
            
            L = 1; % as in Det(0)
        end        
        
        function PH = getPH(self)
            PH = {[-Distrib.InfRate] ,[Distrib.InfRate]};
        end
        
        function bool = isImmediate(self)
            % BOOL = ISIMMEDIATE()
            % Check if the distribution is equivalent to an Immediate
            % distribution
            % Overrides Distrib.isImmediate(self)
            bool = true;
        end        
    end
    
end

