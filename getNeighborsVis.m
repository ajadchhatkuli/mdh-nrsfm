function [ IDX ] = getNeighborsVis( m, Ng, visb )
% get triangulation using pdist functions: Ng number of neighbors
N = length(m(1).m);
distmat = zeros(size(m(1).m,2),size(m(1).m,2),length(m));
for k =1: length(m)
    distmat(:,:,k) = pdist2(m(k).m',m(k).m','cityblock');    
    distmat(~visb(:,k),:,k) = -1;    
%     distmat(:,~visb(:,k),k) = -1;
end

dist = max(distmat,[],3);
% dist = sum(distmat,3);


[~, IDX] = sort(dist,2);

IDX = IDX(:,1:Ng);
% IDX = [(1:N)', IDX];

end