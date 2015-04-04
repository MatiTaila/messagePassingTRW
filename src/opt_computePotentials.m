function [weights, unary, idx] = opt_computePotentials(imL,imR,L,alpha,beta,lmda,K)

PLOT = 0;

imL = double(imL);
imR = double(imR);

[M,N] = size(imL);
nNei  = 4;

% All nodes in image
idx = 1:M*N;

% Take out borders form idx
borders = [1:15*M M+1:M:(N-2)*M+1 (N-1)*M+1:N*M 2*M:M:(N-1)*M];
idx(borders)=[];

% Filter to keep only nodes with label alpha or beta
onlyAlphaBeta = (L(idx)==alpha) | (L(idx)==beta);
idx = idx(onlyAlphaBeta);

% Final quantity of nodes
nNodes = length(idx);

% Find neighbors belonging to alpha or beta
nodes     = repmat(1:nNodes,nNei,1);      nodes     = nodes(:);
neighbors = [idx-1; idx+1; idx-M; idx+M]; neighbors = neighbors(:);
[keep,v2b] = ismember(neighbors,idx);
v1 = nodes(keep);
v2 = v2b(keep);

i=idx(v1);
j=idx(v2);
s = lmda*(1 + (abs(imL(i)-imL(j))<=8));

% Build the weights sparse matrix
weights = sparse(v1,v2,s,nNodes,nNodes); 

% Neighbors that are not in alpha or beta
v1n = nodes(~keep)';
v2n = neighbors(~keep);

tAlpha = zeros(1,nNodes);
tBeta  = zeros(1,nNodes);
for k=1:length(v1n)
    tAlpha(v1n(k))=tAlpha(v1n(k))+min(abs(alpha-L(v2n(k))),K);
    tBeta(v1n(k))=tBeta(v1n(k))+min(abs(beta-L(v2n(k))),K);
end

% Compute unary potentials
alphaInds = idx-alpha*M;
betaInds  = idx-beta*M; 

unary = [ abs(imL(idx)-imR(alphaInds)) + tAlpha;...
    abs(imL(idx)-imR(betaInds)) + tBeta ];

if PLOT
    %%
    figure(100);
    clf
    imagesc(L)
    colormap gray
    hold on
    count = 0;
    for k=1:length(i)
        [p1(2),p1(1)] = ind2sub([M N],i(k));
        [p2(2),p2(1)] = ind2sub([M N],j(k));
        line([p1(1) p2(1)],[p1(2) p2(2)])
        plot(p1(1),p1(2),'*')
        count = count+((L(p1(1),p1(2))==alpha) | (L(p1(1),p1(2))==beta));
    end
%%    
end