classdef Exp < MarkovianDistribution
    % The exponential distribution
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    methods
        function self = Exp(lambda)
            % SELF = EXP(LAMBDA)            
            % Constructs an exponential distribution from the rate
            % parameter
            self@MarkovianDistribution('Exponential', 1);
            setParam(self, 1, 'lambda', lambda, 'java.lang.Double');
            %            self.javaClass = 'jmt.engine.random.Exponential';
            %            self.javaParClass = 'jmt.engine.random.ExponentialPar';
        end
        
        function X = sample(self, n)
            % X = SAMPLE(N)            
            % Get n samples from the distribution
            lambda = self.getParam(1).paramValue;
            X = exprnd(1/lambda,n,1);
        end
        
        function phases = getNumberOfPhases(self)
            % PHASES = GETNUMBEROFPHASES()            
            % Get number of phases in the underpinnning phase-type
            % representation
            phases  = 1;
        end
        
        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)            
            % Evaluate the cumulative distribution function at t
            % AT T
            
            lambda = self.getParam(1).paramValue;
            Ft = 1-exp(-lambda*t);
        end
        
        function PH = getPH(self)
            % PH = GETREPRESENTATION()            
            % Return the renewal process associated to the distribution
            lambda = self.getParam(1).paramValue;
            PH = {[-lambda],[lambda]};
        end
        
        function L = evalLST(self, s)
            % L = EVALST(S)            
            % Evaluate the Laplace-Stieltjes transform of the distribution function at t
            % AT T
            
            lambda = self.getParam(1).paramValue;
            L = lambda / (lambda + s);
        end
        
        function scv = getSCV(self)
            scv = 1.0; 
        end
        
        function update(self,varargin)
            % UPDATE(SELF,VARARGIN)            
            % Update parameters to match the first n central moments
            % (n<=4)
            MEAN = varargin{1};
            SCV = varargin{2};
            SKEW = varargin{3};
            %            KURT = varargin{4};
            if abs(SCV-1) > Distrib.Tol
                line_warning(mfilename,'Warning: the exponential distribution cannot fit squared coefficient of variation != 1, changing squared coefficient of variation to 1.');
            end
            if abs(SKEW-2) > Distrib.Tol
                line_warning(mfilename,'Warning: the exponential distribution cannot fit skewness != 2, changing skewness to 2.');
            end
            %            if abs(KURT-9) > Distrib.Tol
            %                line_warning(mfilename,'Warning: the exponential distribution cannot fit kurtosis != 9, changing kurtosis to 9.');
            %            end
            self.params{1}.paramValue = 1 / MEAN;
        end
        
        function updateMean(self,MEAN)
            % UPDATEMEAN(SELF,MEAN)            
            % Update parameters to match the given mean
            self.params{1}.paramValue = 1 / MEAN;
        end
        
        function updateRate(self,RATE)
            % UPDATERATE(SELF,RATE)            
            % Update rate parameter
            self.params{1}.paramValue = RATE;
        end
        
        function updateMeanAndSCV(self,MEAN,SCV)
            % UPDATEMEANANDSCV(SELF,MEAN,SCV)            
            % Update parameters to match the given mean and squared coefficient of variation (SCV=variance/mean^2)
            if abs(SCV-1) > Distrib.Tol
                line_warning(mfilename,'Warning: the exponential distribution cannot fit SCV != 1, changing SCV to 1.');
            end
            self.params{1}.paramValue = 1 / MEAN;
        end
        
    end
    
    methods (Static)
        function ex = fit(MEAN, SCV, SKEW)
            % EX = FIT(MEAN, SCV, SKEW)
            % Fit the distribution from three standard moments (mean,
            % scv, skewness)
            ex = Exp(1);
            ex.update(MEAN, SCV, SKEW);
        end
        
        function ex = fitMean(MEAN)
            % EX = FITMEAN(MEAN)
            % Fit exponential distribution with given mean
            ex = Exp(1/MEAN);
        end
        
        function ex = fitRate(RATE)
            % EX = FITRATE(RATE)
            % Fit exponential distribution with given rate
            ex = Exp(RATE);
        end
        
        function ex = fitMeanAndSCV(MEAN, SCV)
            % EX = FITMEANANDSCV(MEAN, SCV)
            % Fit exponential distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
            if abs(SCV-1) > Distrib.Tol
                line_warning(mfilename,'Warning: the exponential distribution cannot fit SCV != 1, changing SCV to 1.');
            end
            ex = Exp(1/MEAN);
        end
        
        function Qcell = fromMatrix(Lambda)
            % QCELL = FROMMATRIX(LAMBDA)
            % Instantiates a cell array of Exp objects, each with rate
            % given by the entries of the input matrix
            Qcell = cell(size(Lambda));
            for i=1:size(Lambda,1)
                for j=1:size(Lambda,2)
                    Qcell{i,j} = Exp.fitRate(Lambda(i,j));
                end
            end
        end
    end
    
end