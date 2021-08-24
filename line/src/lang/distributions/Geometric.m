classdef Geometric < DiscreteDistrib
    % A Geometric probability distribution
    %
    % The distribution of the number of Bernoulli trials needed to get 
    % one success.
    %
    % Copyright (c) 2018-2020, Imperial College London
    % All rights reserved.
    
    methods
        function self = Geometric(p)
            % SELF = GEOMETRIC(P)
            self@DiscreteDistrib('Geometric',1,[1,Inf]);
            % Construct a geometric distribution with probability p
            
            setParam(self, 1, 'p', p, 'java.lang.Double');
            %            self.javaClass = '';
            %            self.javaParClass = '';
        end
        
        function ex = getMean(self)
            % EX = GETMEAN()
            
            % Get distribution mean
            p = self.getParam(1).paramValue;

            ex = 1 / p;
        end
        
        function SCV = getSCV(self)
            % SCV = GETSCV()
            
            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            p = self.getParam(1).paramValue;
            
            SCV = 1 - p;
        end
        
        function X = sample(self, n)
            % X = SAMPLE(N)
            if nargin < 2
                n = 1;
            end
            % Get n samples from the distribution
            p = self.getParam(1).paramValue;
            r = rand(n,1);
            X = ceil(log(1-r) ./ log(1-p));
        end
        
        function Ft = evalCDF(self,k)
            % FT = EVALCDF(SELF,K)
            
            % Evaluate the cumulative distribution function at t
            % AT T
            
            p = self.getParam(1).paramValue;
            Ft = 1 - (1-p)^k;
        end
        
        function pr = evalPMF(self, k)
            % PR = EVALPMF(K)
            
            % Evaluate the probability mass function at k
            % AT K
            
            p = self.getParam(1).paramValue;
            pr = (1-p)^(k-1)*p; 
        end
    end
    
end

