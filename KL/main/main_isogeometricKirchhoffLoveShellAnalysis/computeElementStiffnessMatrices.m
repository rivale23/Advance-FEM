function [ BSplinePatch ] = computeElementStiffnessMatrices( BSplinePatch )
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
        
        % check if element is greater than zero
        if Xi(eli+1)~=Xi(eli) && Eta(elj+1)~=Eta(elj)
            %% 3i. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1] 
            %
            %         | xi_i+1 - xi_i                    |
            %         | -------------            0       |
            %         |        2                         |
            %  xi,u = |                                  |
            %         |                  eta_j+1 - eta_j |
            %         |        0         --------------- |
            %         |                          2       |
            detJxiu = (Xi(eli+1)-Xi(eli))*(Eta(elj+1)-Eta(elj))/4;
            
            %% 3ii. Create the Element Freedom Table
            
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
                      
            %% 3iii. Initialize the element area
            elementArea = 0;
            
            %% 3iv. Loop over all the Gauss Points
            
            % initialize element stiffness matrix
            K_local = zeros(noDOFsEl);
            
            for cEta = 1:etaNGP
                for cXi = 1:xiNGP
                    %% 3iv.1. Compute the NURBS coordinates u,v of the Gauss Point coordinates in the bi-unit interval [-1, 1]
                    xi = ( Xi(eli+1)+Xi(eli) + xiGP(cXi)*(Xi(eli+1)-Xi(eli)) )/2;
                    eta = ( Eta(elj+1)+Eta(elj) + etaGP(cEta)*(Eta(elj+1)-Eta(elj)) )/2;
                    
                    %% 3iv.2. Compute the NURBS basis function and up to their second derivatives at the Gauss Point
                    nDrvBasis = 2;
                    dR = computeIGABasisFunctionsAndDerivativesForSurface(eli,p,xi,Xi,elj,q,eta,Eta,CP,isNURBS,nDrvBasis);                    
                    
                    %% 3iv.3. Compute the covariant base vectors and their first derivatives
                    nDrvBaseVct = 1;
                    [dG1,dG2] = computeBaseVectorsAndDerivativesForBSplineSurface(eli,p,elj,q,CP,nDrvBaseVct,dR);
                 
                    %% 3iv.4. Compute the surface normal (third covariant base vector not normalized)
                    G3Tilde = cross(dG1(:,1),dG2(:,1));
                    
                    %% 3iv.5. Compute the legth of G3Tilde (= area dA)
                    dA = norm(G3Tilde);
  
                    %% 3iv.6. Compute the element stiffness matrix at the Gauss point
                    KeOnGP = computeElStiffMtxKirchhoffLoveShellLinear(p,q,dR,[dG1(:,1) dG2(:,1)],[dG1(:,2) dG2(:,2) dG1(:,3)],G3Tilde,Dm,Db);
                    
                    %% 3iv.7 Compute the element area on the Gauss Point and add the contribution
                    elementAreaOnGP = dA*detJxiu*xiGW(cXi)*etaGW(cEta);
                    elementArea = elementArea + elementAreaOnGP;                    
                   
                    %% 3iv.8. Add the contribution from the Gauss Point
                    GaussContribution=KeOnGP*elementAreaOnGP;
                    K_local = K_local + GaussContribution;                    
                end
            end
            %% 3v. Find the minimum element area in the mesh
            if elementArea<minElArea
                minElArea = elementArea;
            end
            
            % save element to BSplinePatch
            BSplinePatch.elements{eli-p,elj-q} = struct('EFT',EFT,...
                                                        'K_local',K_local,...
                                                        'elementArea',elementArea);
        end
    end
end

end

