classdef Transition < Node
    % A service station with queueing
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties
        enablingConditions;
        inhibitingConditions;
        modeNames;
        numbersOfServers;
        timingStrategies;
        distributions;
        firingPriorities;
        firingWeights;
        firingOutcomes;
        % schedStrategy;
        % schedStrategyPar;
        nmodes;
        cap;
        isInit;
    end
    
    methods
        function self = Transition(model,name)
            % TRANSITION(MODEL, NAME)
            
            self@Node(name);
            classes = model.classes;
            self.input = Enabling(classes);
            self.output = Firing(classes);
            self.cap = Inf; % Compatible with other nodes
            
            self.setModel(model);
            self.model.addNode(self);

            self.server = Timing();

            self.enablingConditions = [];
            self.inhibitingConditions = [];
            self.modeNames = {};
            self.numbersOfServers = [];
            self.timingStrategies = [];
            self.distributions = [];
            self.firingPriorities = [];
            self.firingWeights = [];
            self.firingOutcomes = [];
            self.nmodes = 0;
            self.isInit = false;
        end
        
        function bool = isInitialized(self)
            bool = self.isInit;
        end
   
        function self = init(self)
            % SELF = INIT()
            
            nclasses = length(self.model.classes);
            nnodes = length(self.model.nodes);
            self.nmodes = length(self.modeNames);

            self.enablingConditions = zeros(nnodes,nclasses,self.nmodes);
            self.inhibitingConditions = zeros(nnodes,nclasses,self.nmodes);
            self.numbersOfServers = Inf(1,self.nmodes);
            self.timingStrategies = repmat(TimingStrategy.Timed,1,self.nmodes);
            self.firingWeights = ones(1,self.nmodes);
            self.firingPriorities = ones(1,self.nmodes);
            self.distributions = cell(1, self.nmodes);
            self.distributions(:) = {Exp(1)};
            self.firingOutcomes = zeros(nnodes,nclasses,self.nmodes);
            self.isInit = true;
        end
        
        function modeObj = addMode(self, modeName)            
            modeObj = Mode(modeName, self.nmodes+1);
            self.nmodes = self.nmodes + 1;     
            self.modeNames{end+1} = modeName;
            self.enablingConditions(:,end+1) = 0;
            self.inhibitingConditions(:,end+1) = Inf;
            self.numbersOfServers(end+1) = Inf;
            self.timingStrategies(end+1) = TimingStrategy.Timed;
            self.firingWeights(end+1) = 1;
            self.firingPriorities(end+1) = 1;
            self.distributions{end+1} = Exp(1);
            self.firingOutcomes(:,end+1) = 0;
        end

        function self = setEnablingConditions(self, mode, class, node, enablingCondition)
            % SELF = SETENABLINGCONDITIONS(MODE, CLASS, NODE, ENABLINGCONDITIONS)
            
            if ~self.isInitialized
                self.init;
            end
            nnodes = length(self.model.nodes);
            if isa(node, 'Place')
                node = self.model.getNodeIndex(node.name);
                self.enablingConditions(node,class,mode) = enablingCondition;
            elseif isa(node, 'double') && node <= nnodes && isa(self.model.nodes{node}, 'Place')
                self.enablingConditions(node,class,mode) = enablingCondition;
            else
                error('Node must be a Place node or index of a Place node.');
            end
        end

        function self = setInhibitingConditions(self, mode, class, node, inhibitingCondition)
            % SELF = SETINHIBITINGCONDITIONS(MODE, CLASS, NODE, INHIBITINGCONDITIONS)
            
            if ~self.isInitialized
                self.init;
            end
            nnodes = length(self.model.nodes);
            if isa(node, 'Place')
                node = self.model.getNodeIndex(node.name);
                self.inhibitingConditions(node,class,mode) = inhibitingCondition;
            elseif isa(node, 'double') && node <= nnodes && isa(self.model.nodes{node}, 'Place')
                self.inhibitingConditions(node,class,mode) = inhibitingCondition;
            else
                error('Node must be a Place node or index of a Place node.');
            end
        end

        function self = setModeNames(self, mode, modeName)
            % SELF = SETMODENAMES(MODE, MODENAMES)
            
            self.modeNames{mode} = modeName;
        end

        function self = setNumberOfServers(self, mode, numberOfServers)
            % SELF = SETNUMBEROFSERVERS(MODE, NUMOFSERVERS)
            
            self.numbersOfServers(mode) = numberOfServers;
        end

        function self = setTimingStrategy(self, mode, timingStrategy)
            % SELF = SETTIMINGSTRATEGY(MODE, TIMINGSTRATEGY)
            if ~self.isInitialized
                self.init;
            end
            self.timingStrategies(mode) = timingStrategy;
        end
        
        function self = setFiringPriorities(self, mode, firingPriority)
            % SELF = SETFIRINGPRIORITIES(MODE, FIRINGPRIORITIES)
            if ~self.isInitialized
                self.init;
            end            
            self.firingPriorities(mode) = firingPriority;
        end

        function self = setFiringWeights(self, mode, firingWeight)
            % SELF = SETFIRINGWEIGHTS(MODE, FIRINGWEIGHTS)
            if ~self.isInitialized
                self.init;
            end            
            self.firingWeights(mode) = firingWeight;
        end

        function self = setFiringOutcome(self, class, mode, node, firingOutcome)
            % SELF = SETFIRINGOUTCOMES(MODE, CLASS, NODE, FIRINGOUTCOME)
            if ~self.isInitialized
                self.init;
            end
            nnodes = length(self.model.nodes);
            if isa(node, 'Node')
                node = self.model.getNodeIndex(node.name);
                self.firingOutcomes(node,class,mode) = firingOutcome;
            elseif isa(node, 'double') && node <= nnodes && isa(self.model.nodes{node}, 'Node')
                self.firingOutcomes(node,class,mode) = firingOutcome;
            else
                error('Node is not valid.');
            end
        end

        function self = setDistribution(self, mode, distribution)
            self.distributions{mode} = distribution;
        end

        function bool = isTransition(self)
            % BOOL = ISTRANSITION()
            
            bool = isa(self,'Transition');
        end
        
        function sections = getSections(self)
            % SECTIONS = GETSECTIONS()
            
            sections = {self.input, self.server, self.output};
        end


        function [map,mu,phi] = getMarkovianServiceRates(self)
            % [PH,MU,PHI] = GETPHSERVICERATES()
            
            nmodes = self.nmodes;
            map = cell(1,nmodes);
            mu = cell(1,nmodes);
            phi = cell(1,nmodes);

            for r=1:nmodes
                switch class(self.distributions{r})
                    case 'Replayer'
                        aph = self.distributions{r}.fitAPH;
                        map{r} = aph.getRepresentation();
                        mu{r} = aph.getMu;
                        phi{r} = aph.getPhi;
                    case {'Exp','Coxian','Erlang','HyperExp','MarkovianDistribution','APH','MAP'}
                        map{r} = self.distributions{r}.getRepresentation();
                        mu{r} = self.distributions{r}.getMu;
                        phi{r} = self.distributions{r}.getPhi;
                    case 'MMPP2'
                        map{r} = self.distributions{r}.getRepresentation();
                        mu{r} = self.distributions{r}.getMu;
                        phi{r} = self.distributions{r}.getPhi;
                    otherwise
                        map{r}  = {[NaN],[NaN]};
                        mu{r}  = NaN;
                        phi{r}  = NaN;
                end
            end
        end
    end
end
