classdef Host  < LayeredNetworkElement
    % A hardware server in a LayeredNetwork.
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties
        multiplicity;       %int
        replication;        %int
        scheduling;         %char: ps, fcfs, inf, ref
        quantum;            %double
        speedFactor;        %double
        tasks = [];         %list of tasks
    end
    
    methods
        %public methods, including constructor
        
        %constructor
        function self = Host(model, name, multiplicity, scheduling, quantum, speedFactor)
            % self = HOST(MODEL, NAME, MULTIPLICITY, SCHEDULING, QUANTUM, SPEEDFACTOR)
            
            if ~exist('name','var')
                line_error(mfilename,'Constructor requires to specify at least a name.');
            end
            self@LayeredNetworkElement(name);
            
            if ~exist('multiplicity','var')
                multiplicity = 1;
            end
            if ~exist('scheduling','var')
                scheduling = SchedStrategy.PS;
            end
            if ~exist('quantum','var')
                quantum = 0.001;
            end
            if ~exist('speedFactor','var')
                speedFactor = 1;
            end
            self.replication = 1;            
            self.multiplicity = multiplicity;
            self.scheduling = scheduling;
            self.quantum = quantum;
            self.speedFactor = speedFactor;
            model.hosts{end+1} = self;
        end
        
        function self=setReplication(self, replication)
            self.replication = replication;
        end
        
        
        %addTask
        function self = addTask(self, newTask)
            % self = ADDTASK(self, NEWTASK)
            self.tasks = [self.tasks; newTask];            
        end
        
    end
    
end
