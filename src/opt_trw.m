close all
clear all
home

%% Data

imL = convertToGray(imread('../data/tsukuba-imL.png'));
imR = convertToGray(imread('../data/tsukuba-imR.png'));

[M,N] = size(imL);
dMax  = 15;
nLabs = dMax+1;

%% Initializations

% Messages: matrix (nLabs x nNodes). nNodes is different for each tree
mh = zeros(nLabs,N); % horizontal messages
mv = zeros(nLabs,M); % vertical messages

% Unary potentials: array (M x N x nLabs) whith the unary potentials of all
% nodes taking all possible labels. Entry u(i,j,k) is the unary potential
% of pixel (i,j) taking the label k.
u = zeros(M,N,nLabs);

% Pairwise potentials: array (nLabs x nLabs x nNodes). Note that nNodes is
% MxN. The entry pw(i,j,k) is the pairwise potential between nodes k and
% k+1, where node k takes the label i and node k+1 takes the label j. In
% order to be able to define an order k->k+1 we need to keep 2 arrays, one
% for the horizontal passing and one for the vertical.
pwh = zeros(nLabs,nLabs,M*N);
pwv = zeros(nLabs,nLabs,M*N);

% ---------------------- %
% TODO: fill u, pwh, pwv %
% ---------------------- %

%% TRW message passing 

for i=1:M
    for j=1:N
        
        
        
        % unary potentials for node (i,j)
        uij = permute(u(i,j,:),[3 2 1]);
        % update horizontal message
        mh(:,i) = opt_updateMessage(mh(:,i), uij, pwh(:,:,sub2ind([M,N],i,j)));
        % update vertical message
        mv(:,j) = opt_updateMessage(mv(:,j), uij, pwv(:,:,sub2ind([M,N],i,j)));
        
        
        
        
        keyboard
    end
end











% [tree1, tree2] = averageMarginals(tree1,tree2,i,j);

% % Marginals: for each pixel the marginals are a vector in the three
% % dimmensions. i.e.: por pixel (i,j) marginals are m(i,j,:). Marginal
% % asociated with disparity d is m(i,j,d+1) because of lack of index 0
% m = zeros(M,N,nLabs);

% % tree1 and tree2 are the two trees with that variable, with size
% % (dMax+1) x M and (dMax+1) x N respectively
% tree1 = permute(m(i,:,:),[3 2 1]);
% tree2 = permute(m(:,j,:),[3 1 2]);