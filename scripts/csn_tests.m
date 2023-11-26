data = [5 2 3; 3 3 4; 2 3 3; 3 9 1]
[n1,n2] = size(data)
boxsize=0.1
alpha=0.01
weighted=0
c=1:n2

upper = zeros(n1,n2);
lower = zeros(n1,n2);

for i = 1 : n1
% s1 is sorted value of data(i,:);
% s2 is index of element in sort(data(i,:)) corresponding to s1
    [s1,s2] = sort(data(i,:));
% n3 defined as twice of # of cells with negative expression value 
% of specified gene,
% or # of identifier columns.
    n3 = n2-sum(sign(s1));
% h defined as # of cells upper or lower than cell k.
    h = round(boxsize/2*sum(sign(s1)));
    k = 1;
    while k <= n2
% s defined as # of cells with same targeted gene expression as cell k.
        s = 0;
        while k+s+1 <= n2 && s1(k+s+1) == s1(k)
            s = s+1;
        end
% select upper and lower h cells from cell k.
% upper() & lower() are matrices where each element is
% upper or lower threshold of expression level 
% of gene i corresponding to cell k.
        if s >= h
            upper(i,s2(k:k+s)) = data(i,s2(k));
            lower(i,s2(k:k+s)) = data(i,s2(k));
        else
% if k<h | k>n2-h, # of lower/upper cells < h, hence boxsize < 2h.
            upper(i,s2(k:k+s)) = data(i,s2(min(n2,k+s+h)));
% if n3>h & n3+1>k-h, same # of cells as cells with negative values 
% will be eliminated. i.e. Nneg=30, n3-Nneg=60-30=30.
% mistake? 
            lower(i,s2(k:k+s)) = data(i,s2(max(n3*(n3>h)+1,k-h)));
        end
        k = k+s+1;
    end
end

%Construction of cell-specific network
csn = cell(1,n2);
B = zeros(n1,n2);
% idcdf returns inverse cumulative distribution function.
p = -icdf('norm',alpha,0,1);
for k = c
    for j = 1 : n2
        B(:,j) = data(:,j) <= upper(:,k) & data(:,j) >= lower(:,k);
    end
% S = sum(X,DIM) sums along the dimension DIM.
% sum row.
    a = sum(B,2);
    d = (B*B'*n2-a*a')./sqrt((a*a').*((n2-a)*(n2-a)')/(n2-1)+eps);
% interval = n1+1.
    d(1 : n1+1 : end) = 0;
    if weighted
        csn{k} = d.*(d > 0);
    else
% sparse matrix. default.
        csn{k} = sparse(d > p);
    end
    disp(['Cell ' num2str(k) ' is completed']);
end