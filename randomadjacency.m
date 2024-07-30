function [A] = randomadjacency(N,k)
    A = zeros(N,N);
    %Generating Poisson
    outdeg = poissrnd(k,1,N);
    maxdeg = max(outdeg);
    numcombinations = zeros(maxdeg+1,N-1);
    numcombinations(:,1) = 1;
    for i = 2:N-1
        for j = 1:maxdeg+1
            numcombinations(j,i) = sum(numcombinations(1:j,i-1));
        end
    end
    for i = 1:N
        currenttotal = outdeg(i);
        out = zeros(1,N-1);
        for j = 1:N-2
            prob = numcombinations(1:currenttotal+1,N-1-j)'/numcombinations(currenttotal+1,N-j);
            prob = flip(prob);
            cprob = cumsum([0 prob]);
            [~, ~, out(j)] = histcounts(rand, cprob);
            out(j) = out(j)-1;
            currenttotal = currenttotal-out(j);
        end
        out(N-1) = outdeg(i)-sum(out(1:N-2));
        A(:,i) = [out(1:i-1)';0;out(i:N-1)'];
    end
end

