function Offspring = Operator(Parent, SelectedProbability)
% The offspring reproduction operator of 

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    Parent  = Parent.decs;
    Parent1 = Parent(1:floor(end/2), :);
    Parent2 = Parent(floor(end/2)+1:floor(end/2)*2, :);
    N       = size(Parent1, 1);
    
    %% Genetic operators for binary encoding
    Offspring1 = Parent1;
    Offspring2 = Parent2; 
    UnselectedProbability = 1 - SelectedProbability;
    for i = 1 : N
        if rand < 0.5
            index = find(Offspring1(i, :) & ~Offspring2(i, :));
            if size(index, 2) > 0
                % Crossover for offspring1
                k = randperm(size(index, 2), 1);
                index1 = index;
                while k>0
                    ind = index1(TS(UnselectedProbability(index1)));
                    Offspring1(i, ind) = 0;
                    index1(ind==index1) = [];
                    k = k - 1;
                end
                % Crossover for offspring2
                k = randperm(size(index, 2), 1);
                index2 = index;
                while k>0
                    ind = index2(TS(SelectedProbability(index2)));
                    Offspring2(i, ind) = 1;
                    index2(ind==index2) = [];
                    k = k - 1;
                end
            end
            %% Mutation
            % Mutation for offspring1
            index = find(Offspring1(i, :));
            ind = index(TS(UnselectedProbability(index)));
            Offspring1(i, ind) = 0;
            % Mutation for offspring2
            index = find(~Offspring2(i, :));
            ind = index(TS(SelectedProbability(index)));
            Offspring2(i, ind) = 1;
        else
            index = find(~Offspring1(i, :) & Offspring2(i, :));
            if size(index, 2) > 0
                % Crossover for offspring1
                k = randperm(size(index, 2), 1);
                index1 = index;
                while k>0
                    ind = index1(TS(SelectedProbability(index1)));
                    Offspring1(i, ind) = 1;
                    index1(ind==index1) = [];
                    k        = k - 1;
                end
                % Crossover for offspring2
                k = randperm(size(index, 2), 1);
                index2 = index;
                while k>0
                    ind = index2(TS(UnselectedProbability(index2)));
                    Offspring2(i, ind) = 0;
                    index2(ind==index2) = [];
                    k        = k - 1;
                end
            end
            %% Mutation
            % Mutation for offspring1
            index = find(~Offspring1(i, :));
            ind = index(TS(SelectedProbability(index)));
            Offspring1(i, ind) = 1;
            % Mutation for offspring2
            index = find(Offspring2(i, :));
            ind = index(TS(UnselectedProbability(index)));
            Offspring2(i, ind) = 0;
        end
    end
    Offspring = [Offspring1; Offspring2];
    % get unique offspring and individuals (function evaluated)
    Offspring = unique(Offspring, 'rows');
    %% evaluation
    Offspring = SOLUTION(Offspring);
end

function index = TS(Probability)
    % Binary tournament selection based on the probability
    if isempty(Probability)
        index = [];
    else
        index = TournamentSelection(2, 1, -Probability);
    end
end