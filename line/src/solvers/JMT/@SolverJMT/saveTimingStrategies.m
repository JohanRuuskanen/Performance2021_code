function [simDoc, section] = saveTimingStrategies(self, simDoc, section, currentNode)
    % [SIMDOC, SECTION] = SAVETIMINGSTRATEGIES(SIMDOC, SECTION, CURRENTNODE)
    
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.

    strategyNode = simDoc.createElement('parameter');
    strategyNode.setAttribute('array', 'true');
    strategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategy');
    strategyNode.setAttribute('name', 'timingStrategies');
    
    numOfModes = length(currentNode.modeNames);
    for i=1:(numOfModes)
        
        distributionObj = currentNode.distributions{i};
        timingStrategy = currentNode.timingStrategies(i);
        timimgStrategyNode = simDoc.createElement('subParameter');

        if timingStrategy == TimingStrategy.Immediate
            timimgStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy');
            timimgStrategyNode.setAttribute('name', 'DisabledServiceTimeStrategy');
        elseif isa(distributionObj,'APH') || (isa(distributionObj,'Coxian') && distributionObj.getNumParams == 2) || (isa(distributionObj,'HyperExp')  && distributionObj.getNumParams == 2)
            % Coxian and HyperExp have 2 parameters when they have a {mu, p} input specification
            timimgStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy');
            timimgStrategyNode.setAttribute('name', 'timingStrategy');
            distributionNode = simDoc.createElement('subParameter');
            distributionNode.setAttribute('classPath', 'jmt.engine.random.PhaseTypeDistr');
            distributionNode.setAttribute('name', 'Phase-Type');
            distrParNode = simDoc.createElement('subParameter');
            distrParNode.setAttribute('classPath', 'jmt.engine.random.PhaseTypePar');
            distrParNode.setAttribute('name', 'distrPar');
            
            subParNodeAlpha = simDoc.createElement('subParameter');
            subParNodeAlpha.setAttribute('array', 'true');
            subParNodeAlpha.setAttribute('classPath', 'java.lang.Object');
            subParNodeAlpha.setAttribute('name', 'alpha');
            subParNodeAlphaVec = simDoc.createElement('subParameter');
            subParNodeAlphaVec.setAttribute('array', 'true');
            subParNodeAlphaVec.setAttribute('classPath', 'java.lang.Object');
            subParNodeAlphaVec.setAttribute('name', 'vector');
            PH=distributionObj.getRepresentation;
            alpha = abs(map_pie(PH));
            for k=1:distributionObj.getNumberOfPhases
                subParNodeAlphaElem = simDoc.createElement('subParameter');
                subParNodeAlphaElem.setAttribute('classPath', 'java.lang.Double');
                subParNodeAlphaElem.setAttribute('name', 'entry');
                subParValue = simDoc.createElement('value');
                subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',alpha(k))));
                subParNodeAlphaElem.appendChild(subParValue);
                subParNodeAlphaVec.appendChild(subParNodeAlphaElem);
            end
            
            subParNodeT = simDoc.createElement('subParameter');
            subParNodeT.setAttribute('array', 'true');
            subParNodeT.setAttribute('classPath', 'java.lang.Object');
            subParNodeT.setAttribute('name', 'T');
            T = PH{1};
            for k=1:distributionObj.getNumberOfPhases
                subParNodeTvec = simDoc.createElement('subParameter');
                subParNodeTvec.setAttribute('array', 'true');
                subParNodeTvec.setAttribute('classPath', 'java.lang.Object');
                subParNodeTvec.setAttribute('name', 'vector');
                for j=1:distributionObj.getNumberOfPhases
                    subParNodeTElem = simDoc.createElement('subParameter');
                    subParNodeTElem.setAttribute('classPath', 'java.lang.Double');
                    subParNodeTElem.setAttribute('name', 'entry');
                    subParValue = simDoc.createElement('value');
                    if k==j
                        subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',-abs(T(k,j)))));
                    else
                        subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',abs(T(k,j)))));
                    end
                    subParNodeTElem.appendChild(subParValue);
                    subParNodeTvec.appendChild(subParNodeTElem);
                end
                subParNodeT.appendChild(subParNodeTvec);
            end
            
            subParNodeAlpha.appendChild(subParNodeAlphaVec);
            distrParNode.appendChild(subParNodeAlpha);
            distrParNode.appendChild(subParNodeT);
            timimgStrategyNode.appendChild(distributionNode);
            timimgStrategyNode.appendChild(distrParNode);
        elseif (isa(distributionObj,'MAP') && distributionObj.getNumParams == 2)
            timimgStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy');
            timimgStrategyNode.setAttribute('name', 'timingStrategy');
            distributionNode = simDoc.createElement('subParameter');
            distributionNode.setAttribute('classPath', 'jmt.engine.random.MAPDistr');
            distributionNode.setAttribute('name', 'Burst (MAP)');
            distrParNode = simDoc.createElement('subParameter');
            distrParNode.setAttribute('classPath', 'jmt.engine.random.MAPPar');
            distrParNode.setAttribute('name', 'distrPar');
            
            MAP = distributionObj.getRepresentation;
            
            subParNodeD0 = simDoc.createElement('subParameter');
            subParNodeD0.setAttribute('array', 'true');
            subParNodeD0.setAttribute('classPath', 'java.lang.Object');
            subParNodeD0.setAttribute('name', 'D0');
            D0 = MAP{1};
            for k=1:distributionObj.getNumberOfPhases
                subParNodeD0vec = simDoc.createElement('subParameter');
                subParNodeD0vec.setAttribute('array', 'true');
                subParNodeD0vec.setAttribute('classPath', 'java.lang.Object');
                subParNodeD0vec.setAttribute('name', 'vector');
                for j=1:distributionObj.getNumberOfPhases
                    subParNodeD0Elem = simDoc.createElement('subParameter');
                    subParNodeD0Elem.setAttribute('classPath', 'java.lang.Double');
                    subParNodeD0Elem.setAttribute('name', 'entry');
                    subParValue = simDoc.createElement('value');
                    subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',D0(k,j))));
                    subParNodeD0Elem.appendChild(subParValue);
                    subParNodeD0vec.appendChild(subParNodeD0Elem);
                end
                subParNodeD0.appendChild(subParNodeD0vec);
            end
            distrParNode.appendChild(subParNodeD0);
            
            subParNodeD1 = simDoc.createElement('subParameter');
            subParNodeD1.setAttribute('array', 'true');
            subParNodeD1.setAttribute('classPath', 'java.lang.Object');
            subParNodeD1.setAttribute('name', 'D1');
            D1 = MAP{2};
            for k=1:distributionObj.getNumberOfPhases
                subParNodeD1vec = simDoc.createElement('subParameter');
                subParNodeD1vec.setAttribute('array', 'true');
                subParNodeD1vec.setAttribute('classPath', 'java.lang.Object');
                subParNodeD1vec.setAttribute('name', 'vector');
                for j=1:distributionObj.getNumberOfPhases
                    subParNodeD1Elem = simDoc.createElement('subParameter');
                    subParNodeD1Elem.setAttribute('classPath', 'java.lang.Double');
                    subParNodeD1Elem.setAttribute('name', 'entry');
                    subParValue = simDoc.createElement('value');
                    subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',D1(k,j))));
                    subParNodeD1Elem.appendChild(subParValue);
                    subParNodeD1vec.appendChild(subParNodeD1Elem);
                end
                subParNodeD1.appendChild(subParNodeD1vec);
            end
            distrParNode.appendChild(subParNodeD1);
            timimgStrategyNode.appendChild(distributionNode);
            timimgStrategyNode.appendChild(distrParNode);
        else
            timimgStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy');
            timimgStrategyNode.setAttribute('name', 'timingStrategy');
            
            distributionNode = simDoc.createElement('subParameter');
            switch class(distributionObj)
                case 'Det'
                    javaClass = 'jmt.engine.random.DeterministicDistr';
                    javaParClass = 'jmt.engine.random.DeterministicDistrPar';
                case 'Coxian'
                    javaClass = 'jmt.engine.random.CoxianDistr';
                    javaParClass = 'jmt.engine.random.CoxianPar';
                case 'Erlang'
                    javaClass = 'jmt.engine.random.Erlang';
                    javaParClass = 'jmt.engine.random.ErlangPar';
                case 'Exp'
                    javaClass = 'jmt.engine.random.Exponential';
                    javaParClass = 'jmt.engine.random.ExponentialPar';
                case 'Gamma'
                    javaClass = 'jmt.engine.random.GammaDistr';
                    javaParClass = 'jmt.engine.random.GammaDistrPar';
                case 'HyperExp'
                    javaClass = 'jmt.engine.random.HyperExp';
                    javaParClass = 'jmt.engine.random.HyperExpPar';
                case 'Pareto'
                    javaClass = 'jmt.engine.random.Pareto';
                    javaParClass = 'jmt.engine.random.ParetoPar';
                case 'Uniform'
                    javaClass = 'jmt.engine.random.Uniform';
                    javaParClass = 'jmt.engine.random.UniformPar';
                case 'MMPP2'
                    javaClass = 'jmt.engine.random.MMPP2Distr';
                    javaParClass = 'jmt.engine.random.MMPP2Par';
                case 'Replayer'
                    javaClass = 'jmt.engine.random.Replayer';
                    javaParClass = 'jmt.engine.random.ReplayerPar';
            end
            distributionNode.setAttribute('classPath', javaClass);
            switch distributionObj.name
                case 'Replayer'
                    distributionNode.setAttribute('name', 'Replayer');
                case 'HyperExp'
                    distributionNode.setAttribute('name', 'Hyperexponential');
                otherwise
                    distributionNode.setAttribute('name', distributionObj.name);
            end
            timimgStrategyNode.appendChild(distributionNode);
            
            distrParNode = simDoc.createElement('subParameter');
            distrParNode.setAttribute('classPath', javaParClass);
            distrParNode.setAttribute('name', 'distrPar');
            
            for k=1:distributionObj.getNumParams()
                subParNodeAlpha = simDoc.createElement('subParameter');
                subParNodeAlpha.setAttribute('classPath', distributionObj.getParam(k).paramClass);
                subParNodeAlpha.setAttribute('name', distributionObj.getParam(k).paramName);
                subParValue = simDoc.createElement('value');
                switch distributionObj.getParam(k).paramClass
                    case 'java.lang.Double'
                        subParValue.appendChild(simDoc.createTextNode(sprintf('%.12f',distributionObj.getParam(k).paramValue)));
                    case 'java.lang.Long'
                        subParValue.appendChild(simDoc.createTextNode(sprintf('%d',distributionObj.getParam(k).paramValue)));
                    case 'java.lang.String'
                        subParValue.appendChild(simDoc.createTextNode(distributionObj.getParam(k).paramValue));
                end
                subParNodeAlpha.appendChild(subParValue);
                distrParNode.appendChild(subParNodeAlpha);
            end
            
            timimgStrategyNode.appendChild(distrParNode);
        end
        strategyNode.appendChild(timimgStrategyNode);
    end
    section.appendChild(strategyNode);
end
    
    
    