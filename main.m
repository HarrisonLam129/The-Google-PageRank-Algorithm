piprob = [1/4 1/4 1/4 1/4];
d = 0.85;
A = [0 1 0 0;1 0 0 0;1 0 0 1;0 0 0 0];
pi_cprob = cumsum([0, piprob]);
[~, ~, start] = histcounts(rand(1,100),pi_cprob);
%Each column of paths stores a sample path
paths = zeros(1000, 100);
paths(1,:) = start;
for i = 1:1000
    for j = 1:100
        current = paths(i,j);
        if sum(A(:,current)) == 0 %No outgoing links
            [~, ~, paths(i+1,j)] = histcounts(rand, pi_cprob);
        else
            if rand < d
                %Choose an outgoing link at random
                jumpprob = A(:,current)'/sum(A(:,current));
                cprob = cumsum([0, jumpprob]);
                [~, ~, paths(i+1,j)] = histcounts(rand, cprob);
            else
                [~, ~, paths(i+1,j)] = histcounts(rand, pi_cprob);
            end
        end
    end
end

mu = zeros(100,1000,4);
for j = 1:100
    pathj = paths(:,j);
    for k = 1:4
        for t = 1:1000
            mu(j,t,k) = sum(pathj(1:t)==k)/t;
        end
    end
end
hold on
plot(mu(1,:,1))
plot(mu(1,:,2))
plot(mu(1,:,3))
plot(mu(1,:,4))
xlabel('t');
ylabel('$\mu_{1t}^{(k)}$', Interpreter='latex', Rotation=0, Position=[-90 0.5 -1]);
title('$\mu_{1t}^{(k)}$ against t for sample path 1', Interpreter='latex')
legend('k=1', 'k=2', 'k=3', 'k=4')
close

varmu = zeros(4,1000);
for k = 1:4
    for t = 1:1000
        varmu(k,t) = var(mu(:,t,k));
    end
end

hold on
plot(varmu(1,:))
plot(varmu(2,:))
plot(varmu(3,:))
plot(varmu(4,:))
xlabel('t');
ylabel('$Var(\mu_{jt}^{(k)})$', Interpreter='latex', Rotation=0);
title('$Var(\mu_{jt}^{(k)})$ against t', Interpreter='latex')
legend('k=1', 'k=2', 'k=3', 'k=4')
close


scores = pagerank(A,10000);
hold on
plot(mu(1,:,1))
plot(mu(1,:,2))
plot(mu(1,:,3))
plot(mu(1,:,4))
plot(scores(1)/4*ones(1,length(mu(1,:,1))))
plot(scores(2)/4*ones(1,length(mu(1,:,1))))
plot(scores(3)/4*ones(1,length(mu(1,:,1))))
plot(scores(4)/4*ones(1,length(mu(1,:,1))))
xlabel('t');
ylabel('$\mu_{1t}^{(k)}$', Interpreter='latex', Rotation=0, Position=[-90 0.5 -1]);
title('$\mu_{1t}^{(k)}$ against t for sample path 1', Interpreter='latex')
legend('k=1', 'k=2', 'k=3', 'k=4','p1','p2','p3','p4')
close

M = [0 0 5 0;0 0 0 1;5 0 0 0;5 1 5 0];
pagerank(M,10000);


A =  randomadjacency(1000,0.5);
pr = pagerank(A,10000);
pi = ones(1,1000)*1/1000;
newA = zeros(1000);
for i = 1:1000
    if sum(A(:,i)) > 0
        newA(:,i) = A(:,i)*0.85/sum(A(:,i)) + 0.15*pi';
    else
        newA(:,i) = pi';
    end
end
[V,D] = eig(newA);
ev = V(:,1)*1000/sum(V(:,1));
histogram(ev)
title('k = 0.5')
ylabel('Frequency')
xlabel('PageRank score')


links = load('II-9-5-2022-citations.dat');
ids = load('II-9-5-2022-articlejids.dat');
articles1 = unique(links(:));
articles2 = unique(ids(:,1));
%Check if the same set of articles appears in both files
isequal(articles1, articles2);
%Check if articles only appear once in articlejids.dat
isequal(size(ids,1),length(articles2));

z = zeros(1,272);
for i = 1:272
    z(i) = sum(ids(:,2)==i);
end
A = zeros(272);
for i = 1:size(links,1)
    id1 = ids(ids(:,1)==links(i,1),2);
    id2 = ids(ids(:,1)==links(i,2),2);
    A(id2,id1) = A(id2,id1)+1;
end


TC = sum(A,2);
[~,TCRanking] = sort(TC,'descend');
if sum(z) > 0
    pi = z/sum(z);
else
    pi = 0;
end
EF = pagerank(A,10000,pi);
[~,EFRanking] = sort(EF,'descend');
diff1 = zeros(1,272);
for i = 1:272
    diff1(i) = find(TCRanking==i)-find(EFRanking==i);
end

IF = TC./z';
AI = EF./z';
[~,IFRanking] = sort(IF,'descend');
[~,AIRanking] = sort(AI,'descend');
diff2 = zeros(1,272);
for i = 1:272
    diff2(i) = find(IFRanking==i)-find(AIRanking==i);
end





