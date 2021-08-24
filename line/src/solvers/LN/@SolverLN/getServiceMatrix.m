function U = getServiceMatrix(self)
% matrices to obtain entry respt from activities and call respt
eshift = self.lqn.eshift;
U = sparse(self.lqn.nidx + self.lqn.ncalls, self.lqn.nidx + self.lqn.ncalls);
for e = 1:self.lqn.nentries
    eidx = eshift + e;
    U = self.getServiceMatrixRecursion(self.lqn, eidx, U);
end
end

