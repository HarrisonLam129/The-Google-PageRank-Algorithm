function [scores] = pagerank(A,iterations,piprob)
    dim = size(A);
    nodes = dim(1);
    if nargin == 2
        piprob = ones(1,nodes)*(1/nodes);
    end
    d = 0.85;
    S = zeros(nodes);
    %Computing S matrix
    for col = 1:nodes
        if sum(A(:,col)) == 0
            S(:,col) = piprob';
        else
            S(:,col) = d*A(:,col)/sum(A(:,col))+(1-d)*piprob';
        end
    end    
    %Run the recursion
    scores = ones(1,nodes)'*(1/nodes);
    for i = 1:iterations
        scores = S*scores;
    end
    scores = nodes*scores;
end

