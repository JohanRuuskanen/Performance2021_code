classdef Element < Copyable
    % Abstract class for generic elements of a model.
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties
        name;
    end
    
    methods
        %Constructor
        function self = Element(name)
            % SELF = ELEMENT(NAME)
            
            self.setName(name);
        end
        
        function out = getName(self)
            % OUT = GETNAME()
            
            out = self.name;
        end
        
        function self = setName(self, name)
            % SELF = SETNAME(NAME)
            %self.name = categorical({name});
            self.name = name;
        end
        
    end
    
end

