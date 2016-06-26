function [ BSplinePatch] = computeElementStiffnessMatrices( BSplinePatch )
%COMPUTEELEMENTSTIFFNESSMATRICES Summary of this function goes here
%   Detailed explanation goes here

%% 0. Read input

% Assign back the patch properties
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;
parameters = BSplinePatch.parameters;
int = BSplinePatch.int;
TotalArea=0;

% Number of knots in xi-,eta-direction
mxi = length(Xi);
meta = length(Eta);
% Number of control points in xi-,eta-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Create an element freedom table for the patch in the array
BSplinePatch.DOFNumbering = zeros(nxi,neta,3);
k = 1;
for cpj = 1:neta
    for cpi = 1:nxi
        BSplinePatch.DOFNumbering(cpi,cpj,1) = k;
        BSplinePatch.DOFNumbering(cpi,cpj,2) = k + 1;
        BSplinePatch.DOFNumbering(cpi,cpj,3) = k + 2;

        % Update counter
        k = k + 3;
    end
end

% Create the element freedom table for the BSplinePatch into the array of
% the patches
BSplinePatch.EFTPatches = 1:3*BSplinePatch.noCPs;

DOFNumbering = BSplinePatch.DOFNumbering;

% Check input
checkInputForBSplineSurface(p,mxi,nxi,q,meta,neta);
% Initialize minimum element area in the IGA mesh
tolerance = 1e-4;
if abs(CP(1,1,1)-CP(nxi,1,1))>=tolerance
    minElArea = abs(CP(1,1,1)-CP(nxi,1,1));
else
    minElArea = CP(1,1,1)-CP(1,neta,1);
end

% Local number of DOFs per element
noDOFsEl = 3*(p+1)*(q+1);

%% 1. Compute the material matrices

% Compute the membrane material matrix
Dm = parameters.E*parameters.t/(1-parameters.nue^2)*...
    [1              parameters.nue 0
	 parameters.nue 1              0
     0              0              (1-parameters.nue)/2];
                                                 
% Compute the bending material matrix
Db = parameters.E*parameters.t^3/(12*(1-parameters.nue^2))*...
    [1              parameters.nue 0
     parameters.nue 1              0 
     0              0              (1-parameters.nue)/2];

%% 2. Choose an integration rule

% Select the integration scheme
if strcmp(int.type,'default')
    xiNGP = p + 1;
    etaNGP = q + 1;
elseif strcmp(int.type,'manual')
    xiNGP = int.xiNGP;
    etaNGP = int.etaNGP;
end

% Issue the Gauss Point coordinates and weights
[xiGP,xiGW] = getGaussPointsAndWeightsOverUnitDomain(xiNGP);
[etaGP,etaGW] = getGaussPointsAndWeightsOverUnitDomain(etaNGP);

%% 3. loops over elements

n_elements = numel(p+1:mxi-p-1);
m_elements = numel(q+1:meta-q-1);

BSplinePatch.elements = cell(n_elements,m_elements);

for elj = q+1:meta-q-1
    for eli = p+1:mxi-p-1                                    
        %% 3i. Create the Element Freedom Table

        % Initialize element freedom table
        EFT = zeros(1,noDOFsEl);
        % initialize counter
        k = 1;

        % relation global-local dof
        for cpj = elj-q:elj
            for cpi = eli-p:eli
                EFT(k) = DOFNumbering(cpi,cpj,1);
                EFT(k+1) = DOFNumbering(cpi,cpj,2);
                EFT(k+2) = DOFNumbering(cpi,cpj,3);
                % Update counter
                k = k + 3;                                    
            end
        end                                 
        
        %% 3ii. Compute local matrix, element area and total area of the shell
        [ K_local, elementArea ] = computeElementStiffnessMatrix(eli,p,Xi,elj,q,Eta,CP,isNURBS,xiNGP,xiGP,xiGW,etaNGP,etaGP,etaGW,Dm,Db);
        TotalArea=TotalArea+elementArea;

        %% 3iii. save element to BSplinePatch
        element = struct(   'EFT',EFT,...
                            'K_local',K_local,...
                            'elementArea',elementArea);
        BSplinePatch.elements{eli-p,elj-q} = element;  
        
        %% 3iv. Find the minimum element area in the mesh
        if elementArea<minElArea
            minElArea = elementArea;
        end 
    end
end

% save all global quantities to the patch.
BSplinePatch.minElArea = minElArea;
BSplinePatch.TotalArea = TotalArea;
BSplinePatch.TotalMass = TotalArea * parameters.t;

end

