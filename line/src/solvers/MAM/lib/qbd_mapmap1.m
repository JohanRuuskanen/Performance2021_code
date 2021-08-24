function [XN,QN,UN,pqueue,R,eta,G,A_1,A0,A1]=qbd_mapmap1(MAPa,MAPs,util)
% [XN,QN,UN,PQUEUE,R,ETA]=QBD_MAPMAP1(MAPA,MAPS,UTIL)

%[XN,QN,UN,pqueue,R]=qbd_mapmap1(MAPa,MAPs,util)
na = length(MAPa{1});
ns = length(MAPs{1});

if exist('util','var')
    MAPs = map_scale(MAPs,util/map_lambda(MAPa));
end

A1 = kron(MAPa{2},eye(ns));
A0 = krons(MAPa{1},MAPs{1});
A_1 = kron(eye(na),MAPs{2});
A0bar = kron(MAPa{1},eye(ns));

[G,R] = QBD_CR(A_1,A0,A1);

alpha = max(-diag(A0));
A0dt = A0/alpha+eye(size(A0));
A_1dt = A_1/alpha;
A1dt = A1/alpha;
[eta] = QBD_Caudal(A_1dt,A0dt,A1dt);
%%warning off;
pqueue = QBD_pi(A_1,A0bar,R);


if na == 1 && ns == 1
    UN = 1 - pqueue(1);
    QN=(0:(size(pqueue,1)-1))*pqueue;
else
    UN= 1 - sum(pqueue(1,:));
    QN=(0:(size(pqueue,1)-1))*sum(pqueue,2);
end
XN=map_lambda(MAPa);
end
