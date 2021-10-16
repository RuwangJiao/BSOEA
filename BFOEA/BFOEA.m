classdef BFOEA < ALGORITHM
% <multi> <real/binary/permutation> <constrained/none>
% Nondominated sorting genetic algorithm II

%------------------------------- Reference --------------------------------
% K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, A fast and elitist
% multiobjective genetic algorithm: NSGA-II, IEEE Transactions on
% Evolutionary Computation, 2002, 6(2): 182-197.
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
            %% Parameter setting
            LGs = Algorithm.ParameterSet(10);
            
            %% Generate random population
            Population = InitializePopulation(Problem);
            R          = 0.5.*(2.*Problem.D./(2.*Problem.N.*pi)).^(1.0./Problem.D); 
            [~, FrontNo, CrowdDis, ~] = EnvironmentalSelection(Population, Problem.N, R);
            SlidingWindow = [];
            SelectedProbability = zeros(1, Problem.D);
            iter = 0;
            %% Optimization
            while Algorithm.NotTerminated(Population)
                iter = iter + 1;
                radius     = ReduceRadius(R, ceil(Problem.FE/Problem.N), ceil(Problem.maxFE/Problem.N)-1);
                MatingPool = TournamentSelection(2, Problem.N, FrontNo, -CrowdDis);
                Offspring  = OffspringReproduction(Population(MatingPool), SelectedProbability);
                [Population, FrontNo, CrowdDis, NDsolutions] = EnvironmentalSelection([Population, Offspring], Problem.N, radius);            
                [SelectedProbability, SlidingWindow] = CalProbability(iter, LGs, NDsolutions, SlidingWindow);
                
                %%%%% Applied to the test set %%%%%
%                 diversity = DiversityMeasure(Population);
%                 disp(diversity);
                Population = FSTest(Problem, Population);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end
end