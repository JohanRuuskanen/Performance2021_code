classdef JobClassType < NetworkElement
    % An abstract class for a collection of indistinguishable jobs
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        OPEN = 'open';
        CLOSED = 'closed';
        ID_OPEN = 0;
        ID_CLOSED = 1;
    end
    
end
