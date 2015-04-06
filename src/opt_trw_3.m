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
gamma = 0.5;

%% Initializations

% Messages: matrix (nLabs x nNodes). nNodes is different for each tree. The
% entry mh(:,i) corresponds to 
mh = zeros(M,N,nLabs); % Horizontal messages
mv = zeros(M,N,nLabs); % Vertical   messages

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
            
            uHat = uhbar(i,j,:) + mh(i,j,:) + mh(i,j+1,:) + mv(i,j,:) + mv(i+1,j,:);
            
            delta = min(uHat);
            uHat = uHat-delta;
            
            ih = sub2ind([M,N],i,j-1);
            iv = sub2ind([M,N],i-1,j);
            pwhbarij = pwhbar(:,:,ih);
            pwvbarij = pwvbar(:,:,iv);
            
            % Update messages: horizontal
            for k=1:nLabs
                mh(i,j,k) = min( permute(gamma*uHat - mh(i,j-1,:),[3 1 2]) + pwhbarij(:,k) );
            end
            mh(i,j,:) = mh(i,j,k)-min(mh(i,j,k));
            % Update messages: vertical
            for k=1:nLabs
                mv(i,j,k) = min( permute(gamma*uHat - mv(i-1,j,:),[3 1 2]) + pwvbarij(:,k) );
            end
            mv(i,j,:) = mv(i,j,k)-min(mv(i,j,k));
            
        end
    end
    
    for i=M-1:-1:2
        for j=N-1:-1:2
            
            uHat = uhbar(i,j,:) + mh(i,j,:) + mh(i,j+1,:) + mv(i,j,:) + mv(i+1,j,:);
            
            delta = min(uHat);
            uHat = uHat-delta;
            
            ih = sub2ind([M,N],i,j);
            iv = sub2ind([M,N],i,j);
            pwhbarij = pwhbar(:,:,ih)';
            pwvbarij = pwvbar(:,:,iv)';
            
            % Update messages: horizontal
            for k=1:nLabs
                mh(i,j,k) = min( permute(gamma*uHat - mh(i,j,:),[3 1 2]) + pwhbarij(:,k) );
            end
            mh(i,j,:) = mh(i,j,k)-min(mh(i,j,k));
            % Update messages: vertical
            for k=1:nLabs
                mv(i,j,k) = min( permute(gamma*uHat - mv(i,j,:),[3 1 2]) + pwvbarij(:,k) );
            end
            mv(i,j,:) = mv(i,j,k)-min(mv(i,j,k));
            
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
