function jmt
% JMT
% Run java modelling tools start screen

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

if ispc
    cmd = ['java -jar "',jmtGetPath,filesep,'JMT.jar" '];
elseif isunix
    cmd = ['java -jar "',jmtGetPath,filesep,'JMT.jar" '];
else
    cmd = ['java -jar "',jmtGetPath,filesep,'JMT.jar" '];
end
[status] = system(cmd);
if  status > 0
    cmd = ['java --illegal-access=permit -jar "',jmtGetPath,filesep,'JMT.jar" '];
    [status] = system(cmd);
    if status > 0
        rt = java.lang.Runtime.getRuntime();
        rt.exec(cmd);
    end
end
end
