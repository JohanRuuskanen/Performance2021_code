classdef Metric < Copyable
    % An output metric of a Solver, such as a performance index
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        ResidT = 'Residence Time'; % Response Time * Visits
        RespT = 'Response Time'; % Response Time for one Visit
        DropRate = 'Drop Rate';
        QLen = 'Number of Customers';
        QueueT = 'Queue Time';
        FCRWeight = 'FCR Total Weight';
        FCRMemOcc = 'FCR Memory Occupation';
        FJQLen = 'Fork Join Response Time';
        FJRespT = 'Fork Join Response Time';
        RespTSink = 'Response Time per Sink';
        SysDropR = 'System Drop Rate';
        SysQLen = 'System Number of Customers';
        SysPower = 'System Power';
        SysRespT = 'System Response Time';
        SysTput = 'System Throughput';
        Tput = 'Throughput';
        ArvR = 'Arrival Rate';
        TputSink = 'Throughput per Sink';
        Util = 'Utilization';
        TranQLen = 'Tran Number of Customers';
        TranUtil = 'Tran Utilization';
        TranTput = 'Tran Throughput';
        TranRespT = 'Tran Response Time';
        
        ID_ResidT = 0; % Response Time * Visits
        ID_RespT = 1; % Response Time for one Visit
        ID_DropRate = 2;
        ID_QLen = 3;
        ID_QueueT = 4;
        ID_FCRWeight = 5;
        ID_FCRMemOcc = 6;
        ID_FJQLen = 7;
        ID_FJRespT = 8;
        ID_RespTSink = 9;
        ID_SysDropR = 10;
        ID_SysQLen = 11;
        ID_SysPower = 12;
        ID_SysRespT = 13;
        ID_SysTput = 14;
        ID_Tput = 15;
        ID_ArvR = 16;
        ID_TputSink = 17;
        ID_Util = 18;
        ID_TranQLen = 19;
        ID_TranUtil = 20;
        ID_TranTput = 21;
        ID_TranRespT = 22;
    end
    
    
    properties
        type;
        class;
        station;
        simConfInt;
        simMaxRelErr;
        disabled;
        transient;
    end
    
    properties (Hidden)
        stationIndex;
        classIndex;
        %nodeIndex;
    end
    
    methods (Hidden)
        %Constructor
        function self = Metric(type, class, station)
            % SELF = METRIC(TYPE, CLASS, STATION)
            
            self.type = type;
            self.class = class;
            if nargin > 2
                self.station = station;
            else
                self.station = '';
                self.station.name = '';
            end            
            switch self.type
                case {Metric.TranQLen, Metric.TranUtil, Metric.TranTput}
                    self.simConfInt = NaN;
                    self.simMaxRelErr = NaN;
                otherwise % currently used only by JMT
                    self.simConfInt = 0.99;
                    self.simMaxRelErr = 0.03;
            end
            self.disabled = 0;
            self.transient = false;
            switch type
                case {Metric.TranQLen, Metric.TranTput, Metric.TranUtil}
                    self.transient = true;
            end
            self.stationIndex = NaN;
        %    self.nodeIndex = NaN;
            self.classIndex = NaN;
        end
    end
    
    methods
        function self = setTran(self, bool)
            % SELF = SETTRAN(BOOL)
            
            self.transient = bool;
        end
        
        function bool = isTran(self)
            % BOOL = ISTRAN()
            
            bool = self.transient;
        end
        
        function bool = isDisabled(self)
            % BOOL = ISDISABLED()
            
            bool = self.disabled;
        end
        
        function self = disable(self)
            % SELF = DISABLE()
            
            self.disabled = 1;
        end
        
        function self = enable(self)
            % SELF = ENABLE()
            
            self.disabled = 0;
        end
        
        function value = get(self, results, model)
            % VALUE = GET(RESULTS, MODEL)
            
            if self.disabled == 1
                value = NaN;
                return
            end            
            if isnan(self.stationIndex) || self.stationIndex < 0
                stationnames = model.getStationNames();                
                self.stationIndex = findstring(stationnames,self.station.name);
            end
            i = self.stationIndex;
            if isnan(self.classIndex)
                classnames = model.getClassNames();
                self.classIndex = findstring(classnames,self.class.name);
            end
            r = self.classIndex;
            
            switch results.solver
                case 'SolverJMT'                    
                    switch self.type
                        case Metric.TranTput
                            %results.Tran.Avg.T{i,r}.Name = sprintf('Throughput (station %d, class %d)',i,r);
                            %results.Tran.Avg.T{i,r}.TimeInfo.Units = 'since initialization';
                            value = results.Tran.Avg.T{i,r};
                            return
                        case Metric.TranUtil
                            %results.Tran.Avg.U{i,r}.Name = sprintf('Utilization (station %d, class %d)',i,r);
                            %results.Tran.Avg.U{i,r}.TimeInfo.Units = 'since initialization';
                            value = results.Tran.Avg.U{i,r};
                            return
                        case Metric.TranQLen
                            %results.Tran.Avg.Q{i,r}.Name = sprintf('Queue Length (station %d, class %d)',i,r);
                            %results.Tran.Avg.Q{i,r}.TimeInfo.Units = 'since initialization';
                            value = results.Tran.Avg.Q{i,r};
                            return
                        case Metric.TranRespT
                            %results.Tran.Avg.Q{i,r}.Name = sprintf('Queue Length (station %d, class %d)',i,r);
                            %results.Tran.Avg.Q{i,r}.TimeInfo.Units = 'since initialization';
                            value = results.Tran.Avg.R{i,r};
                            return
                    end
                    
                    for i=1:length(results.metric)
                        type = self.type;
                        switch self.type
                            case Metric.TranQLen
                                type = Metric.QLen;
                            case Metric.TranUtil
                                type = Metric.Util;
                            case Metric.TranTput
                                type = Metric.Tput;
                            case Metric.TranRespT
                                type = Metric.RespT;
                        end
                        if strcmp(results.metric{i}.class, self.class.name) && strcmp(results.metric{i}.measureType,type) && strcmp(results.metric{i}.station, self.station.name)
                            chainIdx = find(cellfun(@any,strfind(model.getStruct.classnames,self.class.name)));
                            %chain = model.getChains{chainIdx};
                            switch self.class.type
                                case 'closed'
                                    N = model.getNumberOfJobs();
                                    if results.metric{i}.analyzedSamples > sum(N(chainIdx)) % for a class to be considered recurrent we ask more samples than jobs in the corresponding closed chain
                                        value = results.metric{i}.meanValue;
                                    else
                                        value = 0; % transient metric, long term avg is 0
                                    end
                                case 'open'
                                    if results.metric{i}.analyzedSamples >= 0 % we assume that open classes are always recurrent
                                        value = results.metric{i}.meanValue;
                                    else
                                        value = 0; % transient metric, long term avg is 0
                                    end
                            end
                            break;
                        end
                    end
                otherwise % another LINE solver
                    if ~exist('model','var')
                        line_error(mfilename,'Wrong syntax, use Metric.get(results,model).\n');
                    end
                    if isnan(self.stationIndex)
                        stationnames = model.getStationNames();
                        self.stationIndex = findstring(stationnames,self.station.name);
                    end
                    i = self.stationIndex;
                    if isnan(self.classIndex)
                        classnames = model.getClassNames();
                        self.classIndex = findstring(classnames,self.class.name);
                    end
                    r = self.classIndex;
                    switch self.type
                        case Metric.Util
                            if isempty(results.Avg.U)
                                value = NaN;
                            else
                                value = results.Avg.U(i,r);
                            end
                        case Metric.SysRespT
                            if isempty(results.Avg.C)
                                value = NaN;
                            else
                                value = results.Avg.C(i,r);
                            end
                        case Metric.SysTput
                            if isempty(results.Avg.X)
                                value = NaN;
                            else
                                value = results.Avg.X(i,r);
                            end
                        case Metric.RespT
                            if isempty(results.Avg.R)
                                value = NaN;
                            else
                                value = results.Avg.R(i,r);
                            end                            
                        case Metric.Tput
                            if isempty(results.Avg.T)
                                value = NaN;
                            else
                                value = results.Avg.T(i,r);
                            end
                        case Metric.QLen
                            if isempty(results.Avg.Q)
                                value = NaN;
                            else
                                value = results.Avg.Q(i,r);
                            end
                        case Metric.TranTput
                            %results.Tran.Avg.T{i,r}.Name = sprintf('Throughput (station %d, class %d)',i,r);
                            %results.Tran.Avg.T{i,r}.TimeInfo.Units = 'since initialization';
                            if isempty(results.Tran.Avg.T)
                                value = NaN;
                            else
                                value = results.Tran.Avg.T{i,r};
                            end                            
                        case Metric.TranUtil
                            %results.Tran.Avg.U{i,r}.Name = sprintf('Utilization (station %d, class %d)',i,r);
                            %results.Tran.Avg.U{i,r}.TimeInfo.Units = 'since initialization';
                            value = results.Tran.Avg.U{i,r};
                            if isempty(results.Tran.Avg.U)
                                value = NaN;
                            else
                                value = results.Tran.Avg.U{i,r};
                            end                            
                        case Metric.TranQLen
                            %results.Tran.Avg.Q{i,r}.Name = sprintf('Queue Length (station %d, class %d)',i,r);
                            %results.Tran.Avg.Q{i,r}.TimeInfo.Units = 'since initialization';
                            value = results.Tran.Avg.Q{i,r};
                            if isempty(results.Tran.Avg.Q)
                                value = NaN;
                            else
                                value = results.Tran.Avg.Q{i,r};
                            end                            
                        case Metric.TranRespT
                            %results.Tran.Avg.Q{i,r}.Name = sprintf('Queue Length (station %d, class %d)',i,r);
                            %results.Tran.Avg.Q{i,r}.TimeInfo.Units = 'since initialization';
                            value = results.Tran.Avg.R{i,r};
                            if isempty(results.Tran.Avg.R)
                                value = NaN;
                            else
                                value = results.Tran.Avg.R{i,r};
                            end                            
                    end
            end
        end
    end
end

