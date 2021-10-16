function [Population, FrontNo, CrowdDis, NDsolutions] = EnvironmentalSelection(Population, N, radius)
% The environmental selection of NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    %% calulate niche-count
    niche = zeros(size(Population.objs, 1), 1);
    for i = 1:size(Population.objs, 1)
        for j = i+1:size(Population.objs, 1)
            genes_i = Population(i).decs;
            genes_j = Population(j).decs;
            aDist = pdist2(round(genes_i), round(genes_j), "hamming");
            if aDist < radius
                niche(i) = niche(i) + (1-aDist./radius);
                niche(j) = niche(j) + (1-aDist./radius);
            end
        end
    end
    [FrontNo, MaxFNo] = NDSort([Population.objs, niche], Population.cons, N);
    Next = FrontNo < MaxFNo;
    
    %% Record successful offsprings who are on the first front
    %successful = Population(FrontNo(N+1:end)==1);
    NDsolutions = Population(FrontNo==1);
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population.objs, FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last      = find(FrontNo==MaxFNo);
    [~, Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next); 
end