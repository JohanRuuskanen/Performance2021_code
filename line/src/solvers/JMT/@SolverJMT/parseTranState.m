function [state, evtype, evclass] = parseTranState(fileArv, fileDep, nodePreload)
% [STATE, EVTYPE, EVCLASS] = PARSETRANSTATE(FILEARV, FILEDEP, NODEPRELOAD)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

load(fileArv,'jobArvClasses','jobArvTS','jobArvClassID');
load(fileDep,'jobDepClasses','jobDepTS','jobDepClassID');

%% compute joint state at station
nClasses = length(nodePreload);
state = [jobArvTS,zeros(length(jobArvTS),nClasses);
    jobDepTS,zeros(length(jobDepTS),nClasses)];

evtype = categorical(size(state,1),1);
evclass = zeros(size(state,1),1);

for i=1:size(jobArvTS,1)
    state(i,1+jobArvClassID(i))=+1;
    evtype(i) = EventType.ARV;
    evclass(i) = jobArvClassID(i);
end

for i=1:size(jobDepTS)
    state(length(jobArvTS)+i,1+jobDepClassID(i))=-1;
    evtype(length(jobArvTS)+i) = EventType.DEP;
    evclass(length(jobArvTS)+i) = jobDepClassID(i);
end
[state,I] = sortrows(state,1); % sort on timestamps
state = [0,nodePreload;state];
for j=2:(nClasses+1)
    state(:,j) = cumsum(state(:,j));%+nodePreload(j-1);
end
evtype = [EventType.INIT;evtype(I)']; % categorical requires transpose
evclass = [NaN;evclass(I)];

end