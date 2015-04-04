function m = opt_updateMessage(m,u,pw)
% -------------------------------------------------------------------------
% function m = updateMessage(m,u,pw)
% -------------------------------------------------------------------------
% Computes the update of the message m_{st}. Node s and t are such that s<t
% and the edge s->t is in the set of edges E. In our setting this only
% holds for t=s+1
% -------------------------------------------------------------------------
% inputs:
%   - m: vector (dMax+1) x 1 of messages to be updated, where dMax is the 
%       maximum allowed disparity. Note that (dMax+1=nLabels). This vector
%       encodes the message from s to t: m_{st}
%   - u: vector nLabels x 1 with the unary potentials of node s. Each entry
%       u(i) encodes the unary potential of node s taking the label i.
%   - pw: matrix nLabels x nLabels encoding the pairwise potentials for
%       nodes s and t. Entry pw(i,j) encodes the potential of s taking the
%       label j and t taking the label k.
%   - s,t: nodes involved in the update, such that s < t.
% output:
%   - m: updated message m_{st}
% -------------------------------------------------------------------------

nLabels = size(m,1);
for k=1:nLabels
    m(k) = min( u + pw(:,k) );
end