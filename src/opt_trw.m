close all
clear all
home

%% Data

imL = convertToGray(imread('../data/tsukuba-imL.png'));
imR = convertToGray(imread('../data/tsukuba-imR.png'));

[M,N] = size(imL);
dMax  = 15;
nLabs = dMax+1;
K     = 2;
lmda  = 20;

%% Initializations

% Messages: matrix (nLabs x nNodes). nNodes is different for each tree. The
% entry mh(:,i) corresponds to 
mh = zeros(nLabs,N); % horizontal messages
mv = zeros(nLabs,M); % vertical messages

% Unary potentials: array (M x N x nLabs) whith the unary potentials of all
% nodes taking all possible labels. Entry u(i,j,k) is the unary potential
% of pixel (i,j) taking the label k.

% Pairwise potentials: array (nLabs x nLabs x nNodes). Note that nNodes is
% MxN. The entry pw(i,j,k) is the pairwise potential between nodes k and
% k+1, where node k takes the label i and node k+1 takes the label j. In
% order to be able to define an order k->k+1 we need to keep 2 arrays, one
% for the horizontal passing and one for the vertical.

[uh, uv, pwh, pwv] = opt_computePotentials(imL, imR, nLabs, lmda, K);

%% TRW message passing 

corvengence = 0;
while ~corvengence
   
    for i=2:M-1
        for j=2:N-1
            
            % Unary potentials reparametrization. Here we consider unary
            % potentials for node (i,j)
            uhij = permute(uh(i,j,:),[3 2 1]);
            uvij = permute(uv(i,j,:),[3 2 1]);
            uhij = uhij + mh(:,j);
            uvij = uvij + mv(:,i);
            % Backward message
            uhij = uhij + mh(:,j+1);
            uvij = uvij + mv(:,i+1);


            % Average marginals
            avg = mean([uhij uvij],2);
            uh(i,j,:) = avg;
            uv(i,j,:) = avg;
            
            % Pairwise potentials reparametrization. The -1 is needed in the
            % correct direction (horizontal/vertical) because of the
            % construction of pw. The pairwise between s and t is incoded in
            % position s (here t=s+1, with t being the node (i,j))
            ih = sub2ind([M,N],i,j-1);
            iv = sub2ind([M,N],i-1,j);
            pwhij = pwh(:,:,ih);
            pwvij = pwv(:,:,iv);
            pwhij = bsxfun(@minus, pwhij, mh(:,j));
            pwvij = bsxfun(@minus, pwvij, mv(:,i));
            % Backward message
            pwhij = bsxfun(@minus, pwhij, mh(:,j+1));
            pwvij = bsxfun(@minus, pwvij, mv(:,i+1));
            
            pwh(:,:,ih) = pwhij;
            pwv(:,:,iv) = pwvij;
            
            % update horizontal message
            mh(:,j) = opt_updateMessage(mh(:,j), avg, pwhij);
            % update vertical message
            mv(:,i) = opt_updateMessage(mv(:,i), avg, pwvij);
            
        end
    end
    
    [~,lh] = min(uh,[],3);
    [~,lv] = min(uv,[],3);
    convergence = lh==lv;
    sum(sum(lh~=lv))
    
    figure(1);
    clf;
    imshow(lh,[]);
    drawnow
    
end

% -------------------------------------------------------------------------
% We may need to substract 1 to the final label. disparities goes from 0 to
% 15 and indexes from 1 to 16
% -------------------------------------------------------------------------


figure; 
imshow(lh,[])
