if ~isoctave(), clearvars -except exampleName; end 
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue2', SchedStrategy.PS);

jobclass{1} = ClosedClass(model, 'Class1', 1, node{1}, 0);

servProc1 = Exp(1/0.1);
node{1}.setService(jobclass{1}, servProc1);
servProc2 = Erlang.fitMeanAndSCV(1,1/3);
node{2}.setService(jobclass{1}, servProc2);

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

P = cell(K);
P{1,1} = circul(2);

model.link(P);
%%
solver = SolverFluid(model);
FC = solver.getCdfRespT();
%%
for i=1:model.getNumberOfStations
    for c=1:model.getNumberOfClasses
%        plot(FC{i,c}(:,2),FC{i,c}(:,1)); hold all;
        AvgRespTfromCDF(i,c) = diff(FC{i,c}(:,1))'*FC{i,c}(2:end,2); %mean
        PowerMoment2_R(i,c) = diff(FC{i,c}(:,1))'*(FC{i,c}(2:end,2).^2);
        Variance_R(i,c) = PowerMoment2_R(i,c)-AvgRespTfromCDF(i,c)^2; %variance
        SqCoeffOfVariationRespTfromCDF(i,c) = (Variance_R(i,c))/AvgRespTfromCDF(i,c)^2; %scv
    end
end
fprintf(1,'\n')
disp('Since there is a single job, mean and squared coefficient of variation');
disp('of response times match, up to fluid approximation precision, those');
fprintf(1,'of the service time distribution.\n\n');
AvgRespTfromTheory = [servProc1.getMean; servProc2.getMean]
AvgRespTfromCDF
SqCoeffOfVariationRespTfromTheory = [servProc1.getSCV; servProc2.getSCV]
SqCoeffOfVariationRespTfromCDF
