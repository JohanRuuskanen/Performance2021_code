classdef MMPP2 < MarkovModulated
    % 2-phase Markov-Modulated Poisson Process - MMPP(2)
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    methods
        %Constructor
        function self = MMPP2(lambda0,lambda1,sigma0,sigma1)
            % SELF = MMPP2(LAMBDA0,LAMBDA1,SIGMA0,SIGMA1)
            
            self@MarkovModulated('MMPP2',4);
            setParam(self, 1, 'lambda0', lambda0, 'java.lang.Double');
            setParam(self, 2, 'lambda1', lambda1, 'java.lang.Double');
            setParam(self, 3, 'sigma0', sigma0, 'java.lang.Double');
            setParam(self, 4, 'sigma1', sigma1, 'java.lang.Double');
            %            self.javaClass = 'jmt.engine.random.MMPP2Distr';
            %            self.javaParClass = 'jmt.engine.random.MMPP2Par';
        end
        
        function meant = getMeanT(self,t)
            % MEANT = GETMEANT(SELF,T)
            
            lambda0 =  self.getParam(1).paramValue;
            lambda1 =  self.getParam(2).paramValue;
            sigma0 =  self.getParam(3).paramValue;
            sigma1 =  self.getParam(4).paramValue;
            lambda = (lambda0*sigma1 + lambda1*sigma0) / (sigma0+sigma1);
            meant = lambda * t;
        end
        
        function vart = evalVarT(self,t)
            % VART = EVALVART(SELF,T)
            
            % Evaluate the variance-time curve at t
            lambda0 =  self.getParam(1).paramValue;
            lambda1 =  self.getParam(2).paramValue;
            sigma0 =  self.getParam(3).paramValue;
            sigma1 =  self.getParam(4).paramValue;
            MAP = getRepresentation;
            D0 = MAP{1};
            D1 = MAP{2};
            e = [1;1];
            pie = map_pie(MAP);
            I = eye(2);
            lambda = (lambda0*sigma1 + lambda1*sigma0) / (sigma0+sigma1);
            vart = lambda*t;
            vart = vart + 2*t*(lambda^2-pie*D1*inv(D0+D1+pie*e)*D1*e);
            vart = vart + 2*pi*D1*(expm((D0+D1)*t)-I)*inv(D0+D1+pie*e)^2*D1*e;
        end
        
        % inter-arrival time properties
        function mean = getMean(self)
            % MEAN = GETMEAN()
            
            lambda0 =  self.getParam(1).paramValue;
            lambda1 =  self.getParam(2).paramValue;
            sigma0 =  self.getParam(3).paramValue;
            sigma1 =  self.getParam(4).paramValue;
            lambda = (lambda0*sigma1 + lambda1*sigma0) / (sigma0+sigma1);
            mean = 1 / lambda;
        end
        
        function scv = getSCV(self)
            % SCV = GETSCV()
            
            lambda0 =  self.getParam(1).paramValue;
            lambda1 =  self.getParam(2).paramValue;
            sigma0 =  self.getParam(3).paramValue;
            sigma1 =  self.getParam(4).paramValue;
            scv = (2*lambda0^2*sigma0*sigma1 + lambda0*lambda1*sigma0^2 - 2*lambda0*lambda1*sigma0*sigma1 + lambda0*lambda1*sigma1^2 + lambda0*sigma0^2*sigma1 + 2*lambda0*sigma0*sigma1^2 + lambda0*sigma1^3 + 2*lambda1^2*sigma0*sigma1 + lambda1*sigma0^3 + 2*lambda1*sigma0^2*sigma1 + lambda1*sigma0*sigma1^2)/((sigma0 + sigma1)^2*(lambda0*lambda1 + lambda0*sigma1 + lambda1*sigma0));
        end
        
        function id = getIDC(self)
            % IDC = GETIDC() % INDEX OF DISPERSION FOR COUNTS
            
            lambda0 =  self.getParam(1).paramValue;
            lambda1 =  self.getParam(2).paramValue;
            sigma0 =  self.getParam(3).paramValue;
            sigma1 =  self.getParam(4).paramValue;
            id = 1 + 2*(lambda0-lambda1)^2*sigma0*sigma1/(sigma0+sigma1)^2/(lambda0*sigma1+lambda1*sigma0);
        end
        
        function id = getIDI(self)
            % IDI = GETIDI() % INDEX OF DISPERSION FOR INTERVALS
            
            id = self.getIDC(self); % only asymptotic for now
        end
        
        function MAP = getRepresentation(self)
            MAP = self.getProcess();
        end
        
        function PH = getPH(self)
            % PH = GETREPRESENTATION()
            PH = self.getProcess();
        end
        
        function MAP = getProcess(self)
            % MAP = GETPROCESS()
            
            lambda0 =  self.getParam(1).paramValue;
            lambda1 =  self.getParam(2).paramValue;
            sigma0 =  self.getParam(3).paramValue;
            sigma1 =  self.getParam(4).paramValue;
            D0 = [-sigma0-lambda0,sigma0;sigma1,-sigma1-lambda1];
            D1 = [lambda0,0;0,lambda1];
            MAP = {D0,D1};
        end
        
        function n = getNumberOfPhases(self)
            % N = GETNUMBEROFPHASES()            
            n = 2;
        end
        
        function mu = getMu(self)
            % MU = GETMU()
            MAP = self.getProcess();
            mu = sum(MAP{2},2); % sum D1 rows / diag -D0                        
        end
        
        function phi = getPhi(self)
            % PHI = GETPHI()
            % Return the exit vector of the underlying PH
            MAP = self.getProcess();
            phi = sum(MAP{2},2) ./ diag(-MAP{1}); % sum D1 rows / diag -D0            
        end
        
        function bool = isImmmediate(self)
            % BOOL = ISIMMMEDIATE()
            
            bool = self.getMean() == 0;
        end
    end
end
