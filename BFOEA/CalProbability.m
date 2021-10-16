function [SelectedProbability, SlidingWindow] = CalProbability(generation, LGs, NDsolutions, SlidingWindow)
    %% Calculate the selection frequency probility
    % Remove duplicated solutions in the early generations
    for i = 1:size(SlidingWindow, 2)
        duplicated =  intersect(SlidingWindow(i).ind, NDsolutions.decs, 'row');
        l = SlidingWindow(i).ind;
        for j=1:size(duplicated, 1)
            l(sum(l==duplicated(j,:),2)==size(duplicated, 2),:)= [];
        end
        SlidingWindow(i).ind = l;
    end
    
    % put nondominated solutions into the slidingwindow
    if ceil(generation) > LGs
        for i=1:LGs-1
            SlidingWindow(i).ind = SlidingWindow(i+1).ind;
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
    SelectedProbability = (sum(l,1)./size(l,1));
end