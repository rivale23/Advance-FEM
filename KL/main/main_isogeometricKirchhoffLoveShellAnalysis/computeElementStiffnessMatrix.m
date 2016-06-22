function [ K_local, elementArea ] = computeElementStiffnessMatrix(eli,p,Xi,elj,q,Eta,CP,isNURBS,xiNGP,xiGP,xiGW,etaNGP,etaGP,etaGW,Dm,Db)
%COMPUTEELEMENTSTIFFNESSMATRIX Summary of this function goes here
%   Detailed explanation goes here

% number of DoFs for the element
noDOFsEl = 3*(p+1)*(q+1);
% initialize the element area
elementArea = 0;
% initialize element stiffness matrix
K_local = zeros(noDOFsEl);

% check if element is greater than zero
if Xi(eli+1)~=Xi(eli) && Eta(elj+1)~=Eta(elj)
    %% 1. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1] 
    %
    %         | xi_i+1 - xi_i                    |
    %         | -------------            0       |
    %         |        2                         |
    %  xi,u = |                                  |
    %         |                  eta_j+1 - eta_j |
    %         |        0         --------------- |
    %         |                          2       |
    detJxiu = (Xi(eli+1)-Xi(eli))*(Eta(elj+1)-Eta(elj))/4;         



    %% 2. Loop over all the Gauss Points
    for cEta = 1:etaNGP
        for cXi = 1:xiNGP
            %% 2.1. Compute the NURBS coordinates u,v of the Gauss Point coordinates in the bi-unit interval [-1, 1]
            xi = ( Xi(eli+1)+Xi(eli) + xiGP(cXi)*(Xi(eli+1)-Xi(eli)) )/2;
            eta = ( Eta(elj+1)+Eta(elj) + etaGP(cEta)*(Eta(elj+1)-Eta(elj)) )/2;

            %% 2.2. Compute the NURBS basis function and up to their second derivatives at the Gauss Point
            nDrvBasis = 2;
            dR = computeIGABasisFunctionsAndDerivativesForSurface(eli,p,xi,Xi,elj,q,eta,Eta,CP,isNURBS,nDrvBasis);                    

            %% 2.3. Compute the covariant base vectors and their first derivatives
            nDrvBaseVct = 1;
            [dG1,dG2] = computeBaseVectorsAndDerivativesForBSplineSurface(eli,p,elj,q,CP,nDrvBaseVct,dR);

            %% 2.4. Compute the surface normal (third covariant base vector not normalized)
            G3Tilde = cross(dG1(:,1),dG2(:,1));

            %% 2.5. Compute the legth of G3Tilde (= area dA)
            dA = norm(G3Tilde);

            %% 2.6. Compute the element stiffness matrix at the Gauss point
            KeOnGP = computeElStiffMtxKirchhoffLoveShellLinear(p,q,dR,[dG1(:,1) dG2(:,1)],[dG1(:,2) dG2(:,2) dG1(:,3)],G3Tilde,Dm,Db);

            %% 2.7 Compute the element area on the Gauss Point and add the contribution
            elementAreaOnGP = dA*detJxiu*xiGW(cXi)*etaGW(cEta);
            elementArea = elementArea + elementAreaOnGP;                    

            %% 2.8. Add the contribution from the Gauss Point
            GaussContribution=KeOnGP*elementAreaOnGP;
            K_local = K_local + GaussContribution;                    
        end
    end
end


end

