function ensemble = getEnsemble(self)
% ENSEMBLE = GETENSEMBLE()

if isempty(self.ensemble)
    self.updateEnsemble(true);
end
ensemble = self.ensemble;

end
