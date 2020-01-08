function [mu,D,prob,res,nparams]=NrSfM(IDX,m,vis,solver)
% requires YALMIP functions to call the respective solvers 

% M is the number of views, N the number of tracked points
%%% Inputs:
% IDX is the neighborhood matrix, for each point with K neighbors it is a matrix of
% N x (K+1), on each row the first index is the index of the point and the
% rest the indices of the neighbors

% m is a cell array (1 x M) of normalized point correspondences- each cell is 2 x N.
% vis is a LOGICAL cell array (1 x M) of point visibility, each cell is 1 x N: invisible => 0.

%%% Outputs:
% mu - returned depth matrix at each point on each frame- M x N 
% D  - maximum distance matrix 
% ...
% mu and D are enough 

% solver can be any supported by YALMIP.

% setting default values for visibility (1) and solver (mosek)
if(nargin<4)
    solver='mosek';
    if(nargin<3)
        M=length(m);
        N=size(m{1},2);
        for i=1:M
            vis{i}=ones(N,1);
        end
    end
end


M=length(m);
N=size(m{1},2);
% default parameters for YALMIP SOCP model.
load modelYALMIP;
P=[vis{:}];
P=find(P(:)); % index of visible points


% **: depths for invisible points will be discarded from the
% variable set.
nparams_mu2=length(P); % number of visible points (depths)
P=[P;[N*M+1:N*M+size(IDX,1)*(size(IDX,2)-1)]']; % add index of distance variables (d_ij)
nparams2=length(P); % total number of variables

nparams=N*M+size(IDX,1)*(size(IDX,2)-1); % for size of the constraint matrix (includes invisible 
nparams_mu=N*M; % mu is depth (for indexing purpose)
nparams_D=size(IDX,1)*(size(IDX,2)-1); % D is the distance (for indexing purpose)

IDXt=IDX(:,2:size(IDX,2));  % remove the first point (self)index
IDXt2=repmat(IDX(:,1),1,size(IDXt,2)); % index for the first point in the cone constraint
IDXt3=repmat([1:size(IDX,2)-1],size(IDXt,1),1); % index for the second point in the cone constraint
nconics=0;
for k=1:M
    nconics=nconics+sum(vis{k}(IDXt(:))&vis{k}(IDXt2(:))); % get number of conics
end
A=sparse(nconics*4,nparams); % each cone has four rows, nparams columns and is sparse (it's a big matrix!)
q=4*ones(1,nconics); % store number of rows for each cone
% conic constraints
ni=0;
ni2=0;

for k=1:M % for each view
    pp1=find(vis{k}(IDXt(:))&vis{k}(IDXt2(:))); % get all visible point -neighbor pairs
    nconicsk=numel(pp1); % number of conics for each view
    idxi=IDXt2(pp1); % index for ith (first) point
    idxij=IDXt(pp1); % index for jth (second) point
    idxj=IDXt3(pp1); % index for the d_ij part
    idxk=k*ones(size(idxj)); % view index
    idx1v=sub2ind([N,M],idxi,idxk); % linear index for ith point depth
    idx2v=sub2ind([N,M],idxij,idxk); % linear index for jth point depth
    idx3v=nparams_mu+sub2ind([N,(size(IDX,2)-1)],idxi,idxj); % linear index for d_ij
    allA=[1:4:4*nconicsk]'; % row index for the start of each cone constraint
    qq=[ones(1,nconicsk),m{k}(1,idxi),m{k}(2,idxi),ones(1,nconicsk),-m{k}(1,idxij),-m{k}(2,idxij),-ones(1,nconicsk)]; % **cone constraint
    ind_conic=[allA',allA'+1,allA'+2,allA'+3,allA'+1,allA'+2,allA'+3]; % row index
    ind_param=[idx3v',idx1v',idx1v',idx1v',idx2v',idx2v',idx2v']; % column index
    Ak=sparse(ind_conic,ind_param,qq,nconicsk*4,nparams); % make the constraint matrix for kth view
    A(ni+1:ni+nconicsk*4,:)=Ak;
    ni=ni+nconicsk*4;
    ni2=ni2+nconicsk;
end
A=A(:,P); % remove index of invisible points
%Compose F_struc matrix
C=sparse(nparams2,1); % linear coefficients of the function to be minimized
C(1:nparams_mu2)=-1; % for maximization
A=[sparse(nconics*4,1),[A]]; % add a column for the linear equality constraint (for the constant, see below)
A=[[10,sparse(1,nparams_mu2),-ones(1,nparams_D)];A]; % linear equality constraint sum (D) = 10 (constant)
fprintf('total constraint size = [%d %d]',size(A,1),size(A,2));
model.F_struc=A;
model.c=C;
model.Q=sparse(nparams2,nparams2);
model.monomtable=speye(nparams2);
model.K.f=1; % number of linear equality constraint
model.K.l=0; % number of linear inequalities
model.K.c=0;
model.K.q=q;
model.extended_variables=[];
model.ub=Inf*ones(nparams2,1); % upper bound
model.lb=[sparse(nparams2,1)]; % lower bound
switch(solver)
    case 'mosek'
        output = callmosek(model);
    case 'sedumi'
        output = callsedumi(model);
    otherwise
        output = callmosek(model)
end
solu=output.Primal;
mu=zeros(N,M);
mu(P(1:nparams_mu2))=solu(1:nparams_mu2); % get depth
mu=mu';
D=solu(nparams_mu2+1:nparams2); % get distances
D=reshape(D,size(IDX,1),(size(IDX,2)-1));
end
