function Population = InitializePopulation(Problem)
    Pop = zeros(Problem.N, Problem.D);
    for i = 1 : Problem.N
        k = randperm(round(Problem.D), 1);
        j = randperm(Problem.D, k);
        Pop(i, j) = 1;
    end
    Population = SOLUTION(Pop);
end