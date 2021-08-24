classdef JobClass < NetworkElement
    % An abstract class for a collection of indistinguishable jobs
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties
        priority;
        reference; % reference station
        index;
        type;
        completes; % true if passage through reference station is a completion        
    end
    
    methods (Hidden)
        %Constructor
        function self = JobClass(type, name)
            % SELF = JOBCLASS(TYPE, NAME)
            
            self@NetworkElement(name);
            self.priority = 0;
            self.reference = Node('Unallocated');
            self.index = 1;
            self.type=type;
            self.completes = true;            
        end
        
        function self = setReferenceStation(self, source)
            % SELF = SETREFERENCESTATION(SOURCE)
            
            self.reference = source;
        end
        
        function boolIsa = isReferenceStation(self, node)
            % BOOLISA = ISREFERENCESTATION(NODE)
            
            boolIsa = strcmp(self.reference.name,node.name);
        end
        
        
        %         function self = set.priority(self, priority)
        % SELF = SET.PRIORITY(PRIORITY)
        
        %             if ~(rem(priority,1) == 0 && priority >= 0)
        %                 line_error(mfilename,'Priority must be an integer.\n');
        %             end
        %             self.priority = priority;
        %         end
    end
    
    methods (Access=public)
        function ind = subsindex(self)
            % IND = SUBSINDEX()
            
            ind = double(self.index)-1; % 0 based
        end
    end
    
end
