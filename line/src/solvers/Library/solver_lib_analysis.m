function [QN,UN,RN,TN,CN,XN,runtime] = solver_lib_analysis(qn, options)
% [QN,UN,RN,TN,CN,XN,RUNTIME] = SOLVER_LIB_ANALYSIS(QN, OPTIONS)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

Tstart = tic;

[QN,UN,RN,TN,CN,XN] = solver_lib(qn, options);

QN(isnan(QN))=0;
CN(isnan(CN))=0;
RN(isnan(RN))=0;
UN(isnan(UN))=0;
XN(isnan(XN))=0;
TN(isnan(TN))=0;

runtime = toc(Tstart);

if options.verbose > 0
    line_printf('\nSolver LIB analysis completed. Runtime: %f seconds.\n',runtime);
end
end
