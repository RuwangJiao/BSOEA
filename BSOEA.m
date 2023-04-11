classdef BSOEA < ALGORITHM
% <multi> <binary> <none>
% A Tri-objective Method for Bi-objective Feature Selection in Classification

%------------------------------- Reference --------------------------------
% R. Jiao, B. Xue, and M. Zhang, A Tri-objective method for bi-objective 
% feature selection in classification, Evolutionary Computation (MIT Press), 2023.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting: Learning generations
            LGs = Algorithm.ParameterSet(10);
            
            %% Generate initial population
            Population = InitializePopulation(Problem);
            %% Initialize niche radius
            R                   = 0.5*(2*Problem.D/(2*Problem.N*pi))^(1.0/Problem.D); 
            SlidingWindow       = [];
            SelectedProbability = zeros(1, Problem.D);
            [~, FrontNo, CrowdDis, ~] = EnvironmentalSelection(Population, Problem.N, R);
            iter = 0;
            %% Optimization
            while Algorithm.NotTerminated(Population)
                iter       = iter + 1;
                radius     = ShrinkRadius(R, ceil(Problem.FE/Problem.N), ceil(Problem.maxFE/Problem.N) - 1);
                MatingPool = TournamentSelection(2, Problem.N, FrontNo, -CrowdDis);
                Offspring  = OffspringReproduction(Population(MatingPool), SelectedProbability);
                [Population, FrontNo, CrowdDis, NDsolutions] = EnvironmentalSelection([Population, Offspring], Problem.N, radius);            
                [SelectedProbability, SlidingWindow] = CalProbability(iter, LGs, NDsolutions, SlidingWindow);
            end
        end
    end
end

function [SelectedProbability, SlidingWindow] = CalProbability(generation, LGs, NDsolutions, SlidingWindow)
    %% Calculate the feature selection frequency probility
    % Remove duplicated solutions in the early generations
    for i = 1:size(SlidingWindow, 2)
        duplicated =  intersect(SlidingWindow(i).ind, NDsolutions.decs, 'row');
        l = SlidingWindow(i).ind;
        for j=1:size(duplicated, 1)
            l(sum(l==duplicated(j, :), 2)==size(duplicated, 2), :) = [];
        end
        SlidingWindow(i).ind = l;
    end
    % put nondominated solutions into the slidingwindow
    if ceil(generation) > LGs
        for i = 1:LGs - 1
            SlidingWindow(i).ind = SlidingWindow(i + 1).ind;
        end
        SlidingWindow(LGs).ind = NDsolutions.decs; 
    else
        SlidingWindow(generation).ind = NDsolutions.decs; 
    end
    % Calculate the selection frequency
    l = [];
    for i = 1:size(SlidingWindow, 2)
        l = [l; SlidingWindow(i).ind];
    end
    SelectedProbability = (sum(l, 1)./size(l, 1));
end

function Population = InitializePopulation(Problem)
    %% Initialize an population
    Pop = zeros(Problem.N, Problem.D);
    for i = 1 : Problem.N
        k = randperm(round(Problem.D), 1);
        j = randperm(Problem.D, k);
        Pop(i, j) = 1;
    end
    Population = SOLUTION(Pop);
end

function radius = ShrinkRadius(R, k, MaxK)
    %% The shrink of the dynamic niching radius
    cp       = 5;
    z        = 1e-8;
    Nearzero = 1e-15;
    B        = MaxK./power(log((R + z)./z), 1.0./cp);
    B(B==0)  = B(B==0) + Nearzero;
    f        = R.* exp( -(k./B).^cp );
    tmp      = find(abs(f-z) < Nearzero);
    f(tmp)   = f(tmp).*0 + z;
    radius   = f - z;
    radius(radius<=0) = 0;
end