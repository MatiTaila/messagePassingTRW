function E = opt_computeEnergy(imL,imR,L,K,dMax,lmda)

[M,N] = size(imL);

imL = double(imL);
imR = double(imR);

cols = repmat(dMax+1:N,M,1)-L(:,dMax+1:N);
idxs = (cols-1)*M+repmat((1:M)',1,N-dMax);
unary  = sum(sum(abs(imL(:,dMax+1:N)-imR(idxs))));

lst = (1:M*N)';            % All indexes of matrix in one vector
lst([1:M M+1:M:M*N]) = []; % Take out borders
i = [lst;lst];     i = i(:);
j = [lst-1;lst-M]; j = j(:);
pairwise = sum(lmda*(1+(abs(imL(i)-imL(j))<=8))'*min(abs(L(i)-L(j)),K));

E = unary+pairwise;