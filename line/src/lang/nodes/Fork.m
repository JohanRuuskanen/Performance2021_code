classdef Fork < Node
    % A node to fork jobs into siblings tasks
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties
        schedStrategy;
        cap;
    end
    
    methods
        %Constructor
        function self = Fork(model, name)
            % SELF = FORK(MODEL, NAME)
            
            self@Node(name);
            if(model ~= 0)
                classes = model.classes;
                self.cap = Inf;
                self.input = Buffer(classes);
                self.schedStrategy = SchedStrategy.FORK;
                self.server = ServiceTunnel();
                self.output = Forker(classes);
                self.setModel(model);
                addNode(model, self);
            end
        end
        
        function setTasksPerLink(self, nTasks)
            % SETTASKSPERLINK(NTASKS)
            
            self.output.tasksPerLink = nTasks;
        end
    end
    
end
