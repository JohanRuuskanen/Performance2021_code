classdef ItemSet < NetworkElement
    % A set of cacheable items
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties
        nitems;
        index;
        reference;
        replicable;
    end
    
    methods
        function self = ItemSet(model, name, nitems, reference)
            % SELF = ITEMSET(MODEL, NAME, NITEMS, REFERENCE)
            
            self@NetworkElement(name);
            self.nitems = nitems;
            self.index = 1;
            self.replicable = false;
            if ~isa(reference, 'Cache')
                line_error(mfilename,'ItemClass must be pinned to a Cache.');
            end
            self.reference = reference;
            model.addItemSet(self);
        end
        
        function name = getName(self)
            % NAME = GETNAME()
            
            name = self.name;
        end
        
        function bool = hasReplicableItems(self)
            % BOOL = HASREPLICABLEITEMS()
            
            bool = self.replicable;
        end
        
        function ntypes = getNumberOfItems(self)
            % NTYPES = GETNUMBEROFITEMS()
            
            ntypes = self.nitems;
        end
    end
end
