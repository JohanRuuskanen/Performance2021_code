classdef Place < Station
    % 
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties
        droprule;
        schedStrategies;
        schedStrategy;
        schedStrategyPar;
    end

    methods
        function self = Place(model,name)
            % PLACE(MODEL, NAME)
            
            self@Station(name);
            classes = model.classes;
            self.input = Storage(classes);
            self.output = Linkage(classes);
            self.setModel(model);
            self.model.addNode(self);
            self.server = ServiceTunnel();

            %numOfClasses = [];
            self.numberOfServers = 1;
            self.schedStrategy = SchedStrategy.FCFS;
            self.schedStrategyPar = [];

            self.classCap = [];
            self.cap = [];
            self.droprule = [];
            self.schedStrategies = [];
        end

        function init(self)
            numOfClasses = length(self.model.classes);
            self.schedStrategy = SchedStrategy.FCFS;
            self.schedStrategyPar = zeros(1,numOfClasses);

            self.classCap = Inf(1,numOfClasses);
            self.cap = Inf;
            self.droprule = ones(1, numOfClasses);
            self.schedStrategies = ones(1, numOfClasses);
        end

        function self = setClassCapacity(self, class, capacity)
            % SELF = SETCLASSCAPACITY(CLASS, CAPACITY)
            
            self.classCap(class) = capacity;
        end

        function self = setDropRule(self, class, droprule)
            % SELF = SETDROPRULE(CLASS, DROPRULE)

            self.droprule(class) = droprule;
        end

        function self = setSchedStrategies(self, class, strategy)
            % SELF = SETSCHEDSTRATEGIES(CLASS, STRATEGY)
            
            self.schedStrategies(class) = strategy;
        end

        function sections = getSections(self)
            % SECTIONS = GETSECTIONS()
            
            sections = {self.input, self.server, self.output};
        end
    end
end

