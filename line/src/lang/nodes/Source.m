classdef Source < Station
    % A node to let jobs in open classes arrive to the model
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties
        schedStrategy;
        arrivalProcess;
    end
    
    methods
        %Constructor
        function self = Source(model, name)
            % SELF = SOURCE(MODEL, NAME)
            
            self@Station(name);
            self.numberOfServers = 1;
            if(model ~= 0)
                classes = model.classes;
                self.classCap = Inf*ones(1,length(classes));
                self.output = Dispatcher(classes);
                self.server = ServiceTunnel();
                self.input = RandomSource(classes);
                self.schedStrategy = SchedStrategy.EXT;
                self.setModel(model);
                addNode(model, self);
            end
        end
        
        function self = setScheduling(self, class, strategy)
            %noop
        end
        
        function setArrival(self, class, distribution)
            % SETARRIVAL(CLASS, DISTRIBUTION)
            
            self.input.sourceClasses{1, class.index}{2} = ServiceStrategy.LI;
            self.input.sourceClasses{1, class.index}{3} = distribution;
            self.arrivalProcess{1,class.index} = distribution;
            if distribution.isDisabled()
                self.classCap(class.index) = 0;
            else
                self.classCap(class.index) = Inf;
            end
        end
        
        function sections = getSections(self)
            % SECTIONS = GETSECTIONS()
            
            sections = {self.input, self.server, self.output};
        end
        
        function distrib = getArrivalProcess(self, oclass)
            distrib = self.arrivalProcess{oclass};
        end
        
    end
    
end
