classdef (Sealed) CallType
    properties (Constant)
        SYNC = categorical({'Synchronous'});
        ASYNC = categorical({'Asynchronous'});
        FWD = categorical({'Forwarding'});
        
        ID_SYNC = 1;
        ID_ASYNC = 2;
        ID_FWD = 3;
    end
end