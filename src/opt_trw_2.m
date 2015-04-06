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
gamma = 1;

%% Initializations

% Messages: matrix (nLabs x nNodes). nNodes is different for each tree. The
% entry mh(:,i) corresponds to 
mh = zeros(M,N,nLabs); % horizontal messages
mv = zeros(M,N,nLabs); % vertical messages

% Unary potentials: array (M x N x nLabs) whith the unary potentials of all
% nodes taking all possible labels. Entry u(i,j,k) is the unary potential
% of pixel (i,j) taking the label k.

% Pairwise potentials: array (nLabs x nLabs x nNodes). Note that nNodes is
% MxN. The entry pw(i,j,k) is the pairwise potential between nodes k and
% k+1, where node k takes the label i and node k+1 takes the label j. In
% order to be able to define an order k->k+1 we need to keep 2 arrays, one
% for the horizontal passing and one for the vertical.

[uh, uv, pwh, pwv] = opt_computePotentials(imL, imR, nLabs, lmda, K);
uhbar = uh;
uvbar = uv;
pwhbar = pwh;
pwvbar = pwv;

%% TRW message passing 

corvengence = 0;
while ~corvengence
   
    for i=2:M-1
        for j=2:N-1
            
            % Unary potentials reparametrization. Here we consider unary
            % potentials for node (i,j)
            uhbarij = permute(uhbar(i,j,:),[3 2 1]);
            uvbarij = permute(uvbar(i,j,:),[3 2 1]);
            uhij = uhbarij + permute(mh(i,j,:),[3 2 1]);
            uvij = uvbarij + permute(mv(i,j,:),[3 2 1]);

            % Normalization
            deltah = min(uhij);
            deltav = min(uvij);
            uhij = uhij-deltah;
            uvij = uvij-deltav;
            % Update
            uh(i,j,:) = uhij;
            uv(i,j,:) = uvij;
            
            % Pairwise potentials reparametrization. The -1 is needed in the
            % correct direction (horizontal/vertical) because of the
            % construction of pw. The pairwise between s and t is incoded in
            % position s (here t=s+1, with t being the node (i,j))
            ih = sub2ind([M,N],i,j-1);
            iv = sub2ind([M,N],i-1,j);
%             pwhij = pwh(:,:,ih);
%             pwvij = pwv(:,:,iv);
            pwhbarij = pwhbar(:,:,ih);
            pwvbarij = pwvbar(:,:,iv);
%             pwhij = bsxfun(@minus, pwhij, permute(mh(i,j,:),[3 2 1]));
%             pwvij = bsxfun(@minus, pwvij, permute(mv(i,j,:),[3 2 1]));            
%             % Update
%             pwh(:,:,ih) = pwhij;
%             pwv(:,:,iv) = pwvij;
            
            % Update messages
            mh(i,j,:) = opt_updateMessage(permute(mh(i,j,:),[3 2 1]), uhij, pwhbarij, gamma);
            mh(i,j,:) = opt_updateMessage(permute(mv(i,j,:),[3 2 1]), uvij, pwvbarij, gamma);
            % Normalization
            deltah = min(mh(i,j,:));
            deltav = min(mv(i,j,:));
            mh(i,j,:) = mh(i,j,:)-deltah;
            mv(i,j,:) = mv(i,j,:)-deltav;
            
        end
    end
    
    for i=M-1:-1:2
        for j=N-1:-1:2
            
            % Unary potentials
            uhbarij = permute(uhbar(i,j,:),[3 2 1]);
            uvbarij = permute(uvbar(i,j,:),[3 2 1]);
            uhij = uhbarij + permute(mh(i,j+1,:),[3 2 1]);
            uvij = uvbarij + permute(mv(i+1,j,:),[3 2 1]);
            % Normalization
            deltah = min(uhij);
            deltav = min(uvij);
            uhij = uhij-deltah;
            uvij = uvij-deltav;
            % Update
            uh(i,j,:) = uhij;
            uv(i,j,:) = uvij;
            
            % Pairwise potentials
            ih = sub2ind([M,N],i,j);
            iv = sub2ind([M,N],i,j);
%             pwhij = pwh(:,:,ih);
%             pwvij = pwv(:,:,iv);            
            pwhbarij = pwhbar(:,:,ih);
            pwvbarij = pwvbar(:,:,iv);
%             pwhij = bsxfun(@minus, pwhij, permute(mh(i,j+1,:),[3 2 1]));
%             pwvij = bsxfun(@minus, pwvij, permute(mv(i+1,j,:),[3 2 1]));
%             % Update
%             pwh(:,:,ih) = pwhij;
%             pwv(:,:,iv) = pwvij;
            
            % Update messages
            mh(i,j,:) = opt_updateMessage(permute(mh(i,j,:),[3 2 1]), uhij, pwhbarij', gamma);
            mh(i,j,:) = opt_updateMessage(permute(mv(i,j,:),[3 2 1]), uvij, pwvbarij', gamma);
            % Normalization
            deltah = min(mh(i,j,:));
            deltav = min(mv(i,j,:));
            mh(i,j,:) = mh(i,j,:)-deltah;
            mv(i,j,:) = mv(i,j,:)-deltav;
            
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
