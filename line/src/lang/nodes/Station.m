classdef Station < StatefulNode
    % An abstract class for nodes where jobs station
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties
        numberOfServers;
        cap;
        classCap;
    end
    
    methods(Hidden)
        %Constructor
        function self = Station(name)
            % SELF = STATION(NAME)
            
            self@StatefulNode(name);
            self.cap = Inf;
            self.classCap = [];
        end
    end
    
    methods
        function setNumServers(self, value)
            % SETNUMSERVERS(VALUE)
            
            self.numberOfServers = value;
        end
        
        function setNumberOfServers(self, value)
            % SETNUMBEROFSERVERS(VALUE)
            
            self.numberOfServers = value;
        end
        
        function value = getNumServers(self)
            % VALUE = GETNUMSERVERS()
            
            value = self.numberOfServers;
        end
        
        function value = getNumberOfServers(self)
            % VALUE = GETNUMBEROFSERVERS()
            
            value = self.numberOfServers;
        end
        
        function setCapacity(self, value)
            % SETCAPACITY(VALUE)
            
            self.cap = value;
        end
        
        function setChainCapacity(self, values)
            % SETCHAINCAPACITY(VALUES)
            
            qn = self.model.getStruct;
            if numel(values) ~= qn.nchains
                line_error(mfilename,'The method requires in input a capacity value for each chain.');
            end
            for c = 1:qn.nchains
                inchain = find(qn.chains(c,:));
                for r = inchain
                    if ~self.isServiceDisabled(r)
                        self.classCap(r) = values(c);
                    else
                        self.classCap(r) = Inf;
                    end
                end
            end
            self.cap = min(sum(self.classCap(self.classCap>0)), self.cap);
        end
        
        
        function isD = isServiceDisabled(self, class)
            % ISD = ISSERVICEDISABLED(CLASS)
            
            switch self.server.className
                case 'ServiceTunnel'
                    isD = false;
                otherwise
                    isD = self.server.serviceProcess{1,class}{end}.isDisabled();
            end
        end
        
        function isI = isServiceImmediate(self, class)
            % ISI = ISSERVICEIMMEDIATE(CLASS)
            
            isI = self.server.serviceProcess{1,class}{end}.isImmediate();
        end
        
        function R = getNumberOfServiceClasses(self)
            % R = GETNUMBEROFSERVICECLASSES()
            
            R = size(self.server.serviceProcess,2);
        end
        
        function [p] = getSelfLoopProbabilities(self)
            % [P] = GETSELFLOOPPROBABILITIES()
            
            R = self.getNumberOfServiceClasses();
            p = zeros(1,R);
            for k=1:R
                nOutLinks = length(self.output.outputStrategy{k}{end});
                switch RoutingStrategy.toText(self.output.outputStrategy{k}{2})
                    case 'Random'
                        p(k) = 1 / nOutLinks;
                    case RoutingStrategy.PROB
                        for t=1:nOutLinks % for all outgoing links
                            if strcmp(self.output.outputStrategy{k}{end}{t}{1}.name, self.name)
                                p(k) = self.output.outputStrategy{k}{end}{t}{2};
                                break
                            end
                        end
                end
            end
        end
        
%         function svcProc = getService(self)
%             % svcProc = GETSERVICE()
% 
%             % RETURN SERVICE PROCESSES FOR ALL CLASSES
%             
%             svcProc = {};
%             for r=1:nclasses
%                 if isempty(self.server.serviceProcess{r})
%                     self.server.serviceProcess{r} = {[],ServiceStrategy.LI,Disabled()};
%                     svcProc{r} = Disabled();
%                 elseif self.server.serviceProcess{r}{end}.isImmediate()
%                     svcProc{r} = Immediate();
%                 elseif ~self.server.serviceProcess{r}{end}.isDisabled()
%                     svcProc{r} = serviceProcess{r}{end};
%                 else
%                     svcProc{r} = [];
%                 end
%             end
%         end    

        function [map, mu, phi] = getMarkovianSourceRates(self)
            % [PH,MU,PHI] = GETPHSOURCERATES()
            
            nclasses = size(self.input.sourceClasses,2);
            map = cell(1,nclasses);
            mu = cell(1,nclasses);
            phi = cell(1,nclasses);
            for r=1:nclasses
                if isempty(self.input.sourceClasses{r})
                    self.input.sourceClasses{r} = {[],ServiceStrategy.LI,Disabled()};
                    map{r}  = {[NaN],[NaN]};
                    mu{r}  = NaN;
                    phi{r}  = NaN;
                elseif ~self.input.sourceClasses{r}{end}.isDisabled()
                    switch class(self.input.sourceClasses{r}{end})
                        case 'Replayer'
                            aph = self.input.sourceClasses{r}{end}.fitAPH;
                            map{r} = aph.getRepresentation();
                            mu{r} = aph.getMu;
                            phi{r} = aph.getPhi;
                        case {'Exp','Coxian','Erlang','HyperExp','MarkovianDistribution','APH','PH'}
                            map{r} = self.input.sourceClasses{r}{end}.getRepresentation;
                            mu{r} = self.input.sourceClasses{r}{end}.getMu;
                            phi{r} = self.input.sourceClasses{r}{end}.getPhi;
                        case 'MMPP2'
                            map{r} = self.input.sourceClasses{r}{end}.getRepresentation();
                            mu{r} = self.input.sourceClasses{r}{end}.getMu;
                            phi{r} = self.input.sourceClasses{r}{end}.getPhi;
                        case 'MAP'
                            map{r} = self.input.sourceClasses{r}{end}.getRepresentation();
                            mu{r} = self.input.sourceClasses{r}{end}.getMu;
                            phi{r} = self.input.sourceClasses{r}{end}.getPhi;
                    end
                else
                    map{r}  = {[NaN],[NaN]};
                    mu{r}  = NaN;
                    phi{r}  = NaN;
                end
            end
        end
        
        function [map,mu,phi] = getMarkovianServiceRates(self)
            % [PH,MU,PHI] = GETPHSERVICERATES()
            
            nclasses = size(self.server.serviceProcess,2);
            map = cell(1,nclasses);
            mu = cell(1,nclasses);
            phi = cell(1,nclasses);
            for r=1:nclasses
                if isempty(self.server.serviceProcess{r})
                    self.server.serviceProcess{r} = {[],ServiceStrategy.LI,Disabled()};
                    map{r}  = {[NaN],[NaN]};
                    mu{r}  = NaN;
                    phi{r}  = NaN;
                elseif self.server.serviceProcess{r}{end}.isImmediate()
                    map{r}  = {[-Distrib.InfRate],[Distrib.InfRate]};
                    mu{r}  = [Distrib.InfRate];
                    phi{r}  = [1];
                elseif ~self.server.serviceProcess{r}{end}.isDisabled()
                    switch class(self.server.serviceProcess{r}{end})
                        case 'Replayer'
                            aph = self.server.serviceProcess{r}{end}.fitAPH;
                            map{r} = aph.getRepresentation();
                            mu{r} = aph.getMu;
                            phi{r} = aph.getPhi;
                        case {'Exp','Coxian','Erlang','HyperExp','MarkovianDistribution','APH','MAP'}
                            map{r} = self.server.serviceProcess{r}{end}.getRepresentation();
                            mu{r} = self.server.serviceProcess{r}{end}.getMu;
                            phi{r} = self.server.serviceProcess{r}{end}.getPhi;
                        case 'MMPP2'
                            map{r} = self.server.serviceProcess{r}{end}.getRepresentation();
                            mu{r} = self.server.serviceProcess{r}{end}.getMu;
                            phi{r} = self.server.serviceProcess{r}{end}.getPhi;
                    end
                else
                    map{r}  = {[NaN],[NaN]};
                    mu{r}  = NaN;
                    phi{r}  = NaN;
                end
            end
        end
        
    end
end
