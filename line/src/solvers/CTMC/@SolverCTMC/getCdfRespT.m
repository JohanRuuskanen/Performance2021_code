function RD = getCdfRespT(self, R)
% RD = GETCDFRESPT(R)

if ~exist('R','var')
    R = self.model.getAvgRespTHandles();
end
RD = cell(self.model.getNumberOfStations, self.model.getNumberOfClasses);
qn = self.getStruct;

for i=1:qn.nstations
    for r=1:qn.nclasses
        if ~R{i,r}.disabled
            % tag a class-r job
            taggedModel = self.model.copy;                        
            taggedModel.resetNetwork;
            taggedModel.reset;
            
            taggedModel.classes{end+1,1} = taggedModel.classes{r}.copy;
            if isfinite(qn.njobs(r)) % what if Nr=1 ?
                taggedModel.classes{r}.population = taggedModel.classes{r}.population - 1;
                taggedModel.classes{end,1}.population = 1;
            end
            
            for m=1:length(taggedModel.nodes)
                taggedModel.stations{m}.output.outputStrategy{end+1} = taggedModel.stations{m}.output.outputStrategy{r};
            end
            
            for m=1:length(taggedModel.stations)
                taggedModel.stations{m}.server.serviceProcess{end+1} = taggedModel.stations{m}.server.serviceProcess{r};
                taggedModel.stations{m}.server.serviceProcess{end}{end}=taggedModel.stations{m}.server.serviceProcess{end}{end}.copy;
                if isfinite(taggedModel.stations{m}.classCap)
                    taggedModel.stations{m}.classCap(r) = taggedModel.stations{m}.classCap(r) - 1;
                    taggedModel.stations{m}.classCap(1,end+1) = 1;
                end
            end
            
            %taggedModel.jsimgView
            [Qir,Fir,ev] = SolverCTMC(taggedModel).getGenerator(true);
            
            for v=1:length(ev)
                if ev{v}.passive{1}.event == EventType.ARV && ev{v}.passive{1}.class == r && ev{v}.passive{1}.node == i
                    Air = Fir{v};
                end
            end
            
            for v=1:length(ev)
                if ev{v}.active{1}.event == EventType.DEP && ev{v}.passive{1}.class == r && ev{v}.active{1}.node == i
                    Dir = Fir{v};
                end
            end
            
            A = map_normalize({Qir-Air, Air});
            pie_arv = map_pie(A); % state seen upon arrival of a class-r job
            D = map_normalize({Qir-Dir, Dir});
            
            nonZeroRates = abs(Qir(Qir~=0));
            nonZeroRates = nonZeroRates( nonZeroRates >0 );
            T = abs(100/min(nonZeroRates)); % solve ode until T = 100 events with the slowest rate
            dT = abs(1/max(nonZeroRates)); % solve ode until T = 100 events with the slowest rate
            tset = dT:dT:T;
            
            RD{i,r} = zeros(length(tset),1);
            for t=1:length(tset)
                RD{i,r}(t,2) = tset(t);
                RD{i,r}(t,1) = 1-pie_arv * expm(D{1}*tset(t)) * ones(length(D{1}),1);
            end
        end
    end
end

end