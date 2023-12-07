function [w2m, phim] = FEmodelSolver(N, damagecase, elemass, theta, sen, indexdamage, damageFactor)
% Information for the Model
S = length(sen);         % total number of measured DOFs

% Mass and stiffness Matrix
Mmat    = diag(elemass);

% Formation of stiffness matrix
Kmat       = zeros(N,N); % initial Global stiffness 
elementmatrix   = [1 -1;-1 1]; % element matrix for 2D-problem

% Definition of substructure stiffness matrices
Ksub           = zeros(N,N);
Ksub(:,:,1)    = [1 zeros(1,N-1);zeros(N-1,1) zeros(N-1,N-1)];

for i = 2:N
    Ksub((i-1):i,(i-1):i,i)     = elementmatrix;
end

for j = 1:N
    if (damagecase == 1)&&(indexdamage(j)== 1)
        theta(j) = damageFactor(j)*theta(j);
    end
    Kmat = Kmat + theta(j)*Ksub(:,:,j);
end

% ***************************
%  Generation of modal data
%  ***************************
[modesh0,w2m]   = eig(Kmat,Mmat); % solving the eigen problem
w2m             = diag(w2m); % square frequency or eigen value
[w2m,I]         = sort(w2m);
modesh0         = modesh0(:,I); % modeshape or eigen vector
modesh          = modesh0(sen,:);
phim = zeros(S,N);

for m = 1:N
    modesh(:,m)   = modesh(:,m)/norm(modesh(:,m)); % to normalize all modeshapes to have unit norm
%     modesh(:,m)   = modesh(:,m)/modesh(1,m); % to normalize with respect to top floor level
    phim(:,m) = modesh(:,m)/sign(modesh(end,m));
end