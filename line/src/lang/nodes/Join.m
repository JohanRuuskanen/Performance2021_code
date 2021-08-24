classdef Join < Station
    % A node to join sibling tasks
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    % The number of jobs inside JoinStation is interpreted as the number of
    % jobs waiting to be joined with the sibling tasks.
    
    properties
        joinStrategy;
    end
    
    methods
        %Constructor
        function self = Join(model, name)
            % SELF = JOIN(MODEL, NAME)
            
            self@Station(name);
            if(model ~= 0)
                classes = model.classes;
                self.input = Joiner(classes);
                self.output = Dispatcher(classes);
                self.server = ServiceTunnel();
                self.numberOfServers = Inf;
                self.setModel(model);
                addNode(model, self);
            end
            %             if ~exist('joinstrategy','var')
            %                 joinstrategy = JoinStrategy.STD;
            %             end
            %             setStrategy(joinstrategy);
        end
    end
    
    methods
        function self = setStrategy(self, class, strategy)
            % SELF = SETSTRATEGY(CLASS, STRATEGY)
            
            self.input.setStrategy(class,strategy);
        end
        
        function setScheduling(self, scheduling)
            %noop
        end        
        
        function self = setRequired(self, class, njobs)
            % SELF = SETREQUIRED(CLASS, NJOBS)
            
            self.input.setRequired(class,njobs);
        end
        
        function self = setProbRouting(self, class, destination, probability)
            % SELF = SETPROBROUTING(CLASS, DESTINATION, PROBABILITY)
            
            setRouting(self, class, RoutingStrategy.PROB, destination, probability);
        end
        
    end
    
end
