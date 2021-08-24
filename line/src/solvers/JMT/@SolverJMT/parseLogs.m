function logData = parseLogs(model,isNodeLogged, metric)
% LOGDATA = PARSELOGS(MODEL,ISNODELOGGED, METRIC)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

line_printf('\nJMT log parsing...');
T0=tic;
qn = model.getStruct;
nclasses = qn.nclasses;
logData = cell(qn.nnodes,qn.nclasses);
nodePreload = model.getStateAggr;
for ind=1:qn.nnodes
    if qn.isstateful(ind) && isNodeLogged(ind)
        logFileArv = [model.getLogPath,sprintf('%s-Arv.csv',model.getNodeNames{ind})];
        logFileDep = [model.getLogPath,sprintf('%s-Dep.csv',model.getNodeNames{ind})];
        %% load arrival process
        if exist(logFileArv,'file') && exist(logFileDep,'file')
            %            logArv=readTable(logFileArv,'Delimiter',';','HeaderLines',1); % raw data
            
            %%
            % unclear if this part works fine if user has another local
            % since JMT might write with a different delimiter
            fid=fopen(logFileArv);
            logArv = textscan(fid, '%s%f%f%s%s%s', 'delimiter',';', 'headerlines',1);
            fclose(fid);
            jobArvTS = logArv{2};
            jobArvID = logArv{3};
            jobArvClass = logArv{4};
            clear logArv;
            
            %%
            jobArvClasses = unique(jobArvClass);
            jobArvClassID = zeros(length(jobArvClass),1);
            for c=1:length(jobArvClasses)
                jobArvClassID(find(strcmp(jobArvClasses{c},jobArvClass))) = findstring(model.getClassNames,jobArvClasses{c});
                %                jobArvClassID(find(strcmp(jobArvClasses{c},jobArvClass)))=c;
            end
            logFileArvMat = [model.getLogPath,filesep,sprintf('%s-Arv.mat',model.getNodeNames{ind})];
            save(logFileArvMat,'jobArvTS','jobArvID','jobArvClass','jobArvClasses','jobArvClassID');
            
            %% load departure process
            fid=fopen(logFileDep);
            logDep = textscan(fid, '%s%f%f%s%s%s', 'delimiter',';', 'headerlines',1);
            fclose(fid);
            jobDepTS = logDep{2};
            jobDepID = logDep{3};
            jobDepClass = logDep{4};
            clear logDep;
            
            %             jobDepTS = table2array(logDep(:,2));
            %             jobDepID = table2array(logDep(:,3));
            %             jobDepClass = table2cell(logDep(:,4));
            jobDepClasses = unique(jobDepClass);
            jobDepClassID = zeros(length(jobDepClass),1);
            for c=1:length(jobDepClasses)
                jobDepClassID(find(strcmp(jobDepClasses{c},jobDepClass))) = findstring(model.getClassNames,jobDepClasses{c});
                %                jobDepClassID(find(strcmp(jobDepClasses{c},jobDepClass)))=c;
            end
            logFileDepMat = [model.getLogPath,filesep,sprintf('%s-Dep.mat',model.getNodeNames{ind})];
            save(logFileDepMat,'jobDepTS','jobDepID','jobDepClass','jobDepClasses','jobDepClassID');
            
            switch metric
                case Metric.QLen
                    [nodeState{ind}, evtype, evclass] = SolverJMT.parseTranState(logFileArvMat, logFileDepMat, nodePreload{ind});
                    
                    %% save in default data structure
                    for r=1:nclasses %0:numOfClasses
                        logData{ind,r} = struct();
                        logData{ind,r}.t = nodeState{ind}(:,1);
                        logData{ind,r}.event = cell(length(evtype),1);
                        ec = 0;
                        for e=1:length(evtype)
                            if evclass(e) == r
                                logData{ind,r}.event{e,1} = Event(evtype(e), ind, r, NaN, [], nodeState{ind}(e,1));
                            else
                                logData{ind,r}.event{e,1} = [];
                            end
                        end
                        logData{ind,r}.QLen = nodeState{ind}(:,1+r);
                        logData{ind,r}.arvID = jobArvID;
                        logData{ind,r}.depID = jobDepID;                        
                    end
                case Metric.RespT
                    [classResT, jobRespT, jobResTArvTS] = SolverJMT.parseTranRespT(logFileArvMat, logFileDepMat);
                    
                    for r=1:nclasses
                        logData{ind,r} = struct();
                        if r <= size(classResT,2)
                            logData{ind,r}.t = jobResTArvTS;
                            logData{ind,r}.RespT = classResT{r};
                            %logData{i,r}.PassT = jobRespT;
                        else
                            logData{ind,r}.t = [];
                            logData{ind,r}.RespT = [];
                            %logData{i,r}.PassT = [];
                        end
                    end
            end
        end
    end
end
runtime=toc(T0);
line_printf(' completed. Runtime: %f seconds.\n',runtime);
end

