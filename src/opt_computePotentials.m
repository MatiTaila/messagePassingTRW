function [uh, uv, pwh, pwv] = opt_computePotentials(imL, imR, nLabs, lmda, K)
% -------------------------------------------------------------------------
% function [uh, uv, pwh, pwv] = opt_computePotentials(imL, imR, nLabs,
%   lmda, K);
% -------------------------------------------------------------------------
% Unary potentials: array (M x N x nLabs) whith the unary potentials of all
% nodes taking all possible labels. Entry u(i,j,k) is the unary potential
% of pixel (i,j) taking the label k.
% uh = zeros(M,N,nLabs);
% uv = zeros(M,N,nLabs);
% -------------------------------------------------------------------------
% Pairwise potentials: array (nLabs x nLabs x nNodes). Note that nNodes is
% MxN. The entry pw(i,j,k) is the pairwise potential between nodes k and
% k+1, where node k takes the label i and node k+1 takes the label j. In
% order to be able to define an order k->k+1 we need to keep 2 arrays, one
% for the horizontal passing and one for the vertical.
% pwh = zeros(nLabs,nLabs,M*N);
% pwv = zeros(nLabs,nLabs,M*N);
% -------------------------------------------------------------------------

[M,N] = size(imL);

% Pairwise potentials
Vpq = K*ones(nLabs);
Vpq = Vpq-diag(K*ones(nLabs,1));
for i=1:K-1
    Vpq = Vpq-diag((K-i)*ones(nLabs-i,1),i);
    Vpq = Vpq-diag((K-i)*ones(nLabs-i,1),-i);
end
Vpqh = repmat(Vpq,[1 1 M*N]);
Vpqv = repmat(Vpq,[1 1 M*N]);

wpqh = 1 + (abs([imL(:,2:N)-imL(:,1:N-1) zeros(M,1)])<=8);
wpqv = 1 + (abs([imL(2:M,:)-imL(1:M-1,:);zeros(1,N)])<=8);

pwh = lmda*permute(repmat(wpqh(:),[1,nLabs,nLabs]),[2 3 1]).*Vpqh;
pwv = lmda*permute(repmat(wpqv(:),[1,nLabs,nLabs]),[2 3 1]).*Vpqv;

% Unary potentials
uh = zeros(M,N,nLabs);
uv = zeros(M,N,nLabs);
for i=1:nLabs
    d = i-1;
    uh(:,:,i) = abs(imL - [zeros(M,d) imR(:,1:end-d)]);
    uv(:,:,i) = abs(imL - [zeros(d,N);imR(1:end-d,:)]);
end
