function [ K_global, F_global ] = assembleGlobalSystem( BSplinePatch )
%ASSEMBLEGLOBALSYSTEM Summary of this function goes here
%   Detailed explanation goes here

%% 0. Read input

% Assign back the patch properties
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;
int = BSplinePatch.int;

% Neuman boundary conditions
NBC = BSplinePatch.NBC;

% Number of knots in xi-,eta-direction
mxi = length(Xi);
meta = length(Eta);
% Number of control points in xi-,eta-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Check input
checkInputForBSplineSurface(p,mxi,nxi,q,meta,neta);

%% 1. loops over elements

% Global number of DOFs
noDOFs = 3*nxi*neta;
K_global = zeros(noDOFs);

for elj = q+1:meta-q-1
    for eli = p+1:mxi-p-1
        % check if element is greater than zero
        if Xi(eli+1)~=Xi(eli) && Eta(elj+1)~=Eta(elj)
            
            element = BSplinePatch.elements{eli-p,elj-q};
            
            %% 1i. Read the Element Freedom Table            
            EFT = element.EFT;
            
            %% 1ii. Add element contribution to global stiffness matrix
            
            K_local = element.K_local;
            K_global(EFT,EFT) = K_global(EFT,EFT) + K_local;
            
        end
    end
end

%% 2. Compute the exernally applied load vector
F_global = zeros(noDOFs,1);
t = 0;    
for counterNBC = 1:NBC.noCnd
    funcHandle = str2func(NBC.computeLoadVct{counterNBC});
    F_global = funcHandle(F_global,NBC.xiLoadExtension{counterNBC},...
        NBC.etaLoadExtension{counterNBC},p,q,Xi,Eta,CP,isNURBS,...
        NBC.loadAmplitude{counterNBC},...
        NBC.loadDirection(counterNBC,1),t,int,'');
end

%% 3. Update the right-hand side vector if inhomogeneous Dirichlet boundary conditions are encountered

% Get the numbering and the values of the DOFs which are prescribed
inhomDOFs = BSplinePatch.inhomDOFs;
valuesInhomDOFs = BSplinePatch.valuesInhomDOFs;

if norm(valuesInhomDOFs(inhomDOFs)) ~= 0
    F_global = F_global - K_global(:,inhomDOFs)*valuesInhomDOFs(inhomDOFs);
end

end

