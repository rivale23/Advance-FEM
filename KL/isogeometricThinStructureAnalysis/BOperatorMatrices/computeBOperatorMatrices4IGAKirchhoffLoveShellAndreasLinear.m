%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BOperatorMembraneStrainContra,BOperatorBendingStrainContra,...
    dOperatorBendingStraindxiContra,dOperatorBendingStraindetaContra] = ...
    computeBOperatorMatrices4IGAKirchhoffLoveShellAndreasLinear...
    (dRdxi,dRdeta,ddRdxidxi,ddRdetadeta,ddRdxideta,dddRdxidxidxi,...
    dddRdxidxideta,dddRdetadetadxi,dddRdetadetadeta,dACovariant1,...
    dACovariant2,A3Tilde,dA3dxi,dA3deta,dAKommaAlpha,BV)
%% Function documentation
%
% Returns the B-operator matrices for the membrane and the bending strain
% Voigt vectors as well as for the derivatives of the bending strain with 
% respect to both parametric directions in the contravariant basis. That 
% is, multiplication of the those matrices with the vector of DOFs returns 
% the membrane and the bending strain vectors as well as the derivatives of 
% the bending strain vector with respect to both parametric directions.
%
%                            Input :
%                              p,q : polynomial degrees
%                            dRdxi : Matrix containing the derivatives of
%                                    the basis functions dR/dxi
%                           dRdeta : Matrix containing the derivatives of
%                                    the basis functions dR/deta
%                        ddRdxidxi : Matrix containing the derivatives of
%                                    the basis functions d^2R/dxi^2
%                      ddRdetadeta : Matrix containing the derivatives of
%                                    the basis functions d^2R/deta^2
%                       ddRdxideta : Matrix containing the derivatives of
%                                    the basis functions d^2R/dxi/deta
%                    dddRdxidxidxi : Matrix containing the derivatives of
%                                    the basis functions d^3R/dxi^3
%                   dddRdxidxideta : Matrix containing the derivatives of
%                                    the basis functions d^3R/dxi^2/deta
%                 dddRdetadetadeta : Matrix containing the derivatives of
%                                    the basis functions d^3R/deta^3
%                     dACovariant1 : The base vector A1 and up to its 
%                                    second mixed derivatives :
%                                    = [A1 dA1/dxi d^2A1/dxi^2 
%                                       dA1/deta(=dA2/dxi) 
%                                       d^2A1/detadxi(=d^2A2/dxi^2) 
%                                       d^2A1/deta^2(=d^2A2/dxideta)]
%                     dACovariant2 : The base vector A1 and up to its 
%                                    second mixed derivatives (symmetry 
%                                    with respect to the A1 base vector is 
%                                    taken into consideration)
%                                    = [A2 dA2/deta d^2A2/deta^2]
%                          A3Tilde : The not normalized surface normal
%                           dA3dxi : The parametric derivative of the 
%                                    surface normal dA3/dxi
%                          dA3deta : The parametric derivative of the 
%                                    surface normal dA3/deta
%                     dAKommaAlpha : The derivatives of the element surface
%                                    area
%                               BV : The components of the curvature tensor
%                                    in Voigt notation = [B11 B22 B12]
%
%                           Output :
%    BOperatorMembraneStrainContra : The B-operator matrix for the membrane 
%                                    strain in the contravariant basis 
%     BOperatorBendingStrainContra : The B-operator matrix for the bending 
%                                    strain in the contravariant basis
%  dOperatorBendingStraindxiContra : The B-operator for the computation of
%                                    the derivative in the xi-parametric
%                                    direction of the bending strain
% dOperatorBendingStraindetaContra : The B-operator for the computation of
%                                    the derivative in the eta-parametric
%                                    direction of the bending strain
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the parametric derivatives of the not normalized surface normal
%
% 2. Compute the D matrix (see Nitsche - Computation Vol. 4) which is the variation of W vector (see equation (1c) Nitsche - Computation Vol. 5) 
%
% 3. Compute the variation of the C vector (see Nitsche - Computation Vol. 4)
%
% 4. Compute variations of vectors W,gamma (see equation (8) Nitsche - Computation Vol. 5)
%
% 5. Compute the variations of the derivatives of the strain components in the contravariant basis (see equation (2) Nitsche - Computation Vol. 5)
%
% 6. Compute the B-operator matrices for the membrane strain, the bending strain as well as the parametric derivatives of the bending strain in the covariant basis
%
%% Function main body

%% 0. Read input

% For the permutation matrices needed for the cross products

% A1Cross Matrix
A1Cross = [0                   -dACovariant1(3,1) dACovariant1(2,1)
           dACovariant1(3,1)   0                  -dACovariant1(1,1)
           -dACovariant1(2,1)  dACovariant1(1,1)  0                 ];
   
% A2Cross Matrix
A2Cross = [0                   -dACovariant2(3,1) dACovariant2(2,1)
           dACovariant2(3,1)   0                  -dACovariant2(1,1)
           -dACovariant2(2,1)  dACovariant2(1,1)  0                 ];
   
% dA2dxiCross Matrix
dA2dxiCross = [0                   -dACovariant1(3,4) dACovariant1(2,4)
               dACovariant1(3,4)   0                  -dACovariant1(1,4)
               -dACovariant1(2,4)  dACovariant1(1,4)  0                 ];
           
% dA1dxiCross Matrix 
dA1dxiCross = [0                   -dACovariant1(3,2) dACovariant1(2,2)
               dACovariant1(3,2)   0                  -dACovariant1(1,2)
               -dACovariant1(2,2)  dACovariant1(1,2)  0                 ];

% dA2detaCross Matrix
dA2detaCross = [0                   -dACovariant2(3,2) dACovariant2(2,2)
                dACovariant2(3,2)   0                  -dACovariant2(1,2)
                -dACovariant2(2,2)  dACovariant2(1,2)  0                 ];
            
% Surface normal base vector
A3 = A3Tilde/norm(A3Tilde);

%% Compute the parametric derivatives of the curvature tensor

% dBV11dalpha
dBV11dxi = - dACovariant1(:,3)'*A3 - dACovariant1(:,2)'*dA3dxi;
dBV11deta = - dACovariant1(:,5)'*A3 - dACovariant1(:,2)'*dA3deta;

% dBV22dalpha
dBV22dxi = - dACovariant1(:,6)'*A3 - dACovariant2(:,2)'*dA3dxi;
dBV22deta = - dACovariant2(:,3)'*A3 - dACovariant2(:,2)'*dA3deta;

% dBV12dalpha
dBV12dxi = - dACovariant1(:,5)'*A3 - dACovariant1(:,4)'*dA3dxi;
dBV12deta = - dACovariant1(:,6)'*A3 - dACovariant1(:,4)'*dA3deta;

% Compute the variation of the C scalar as it is defined in equation (5c) -
% Nitsche Computation Vol. 10
CKommaU = cross(dACovariant2(:,1),A3)'*dRdxi - cross(dACovariant1(:,1),A3)'*dRdeta;

%% Compute the inverse of the element area
jInv = 1/norm(A3Tilde);

%% Compute the parametric derivatives of the inverse of the element surface area, see equation (2) - Nitsche Computation Vol. 10

% dj^{-1}dxi
djInvdxi = -1/(norm(A3Tilde))*dAKommaAlpha(1,1);

% dj^{-1}dxi
djInvdeta = -1/(norm(A3Tilde))*dAKommaAlpha(2,1);

%% Compute the variations of the scalar C as in equation (5d) - Nitsche  Computation Vol. 10

% dC/dxi
CKomma1KommaU = (cross(dACovariant1(:,4),A3) + cross(dACovariant2(:,1),dA3dxi))'*dRdxi + ...
    cross(dACovariant2(:,1),A3)'*ddRdxidxi - ...
    (cross(dACovariant1(:,2),A3) + cross(dACovariant1(:,1),dA3dxi))'*dRdeta - ...
    cross(dACovariant1(:,1),A3)'*ddRdxideta;

% dC/deta
CKomma2KommaU = (cross(dACovariant2(:,2),A3) + cross(dACovariant2(:,1),dA3deta))'*dRdxi + ...
    cross(dACovariant2(:,1),A3)'*ddRdxideta - ...
    (cross(dACovariant1(:,4),A3) + cross(dACovariant1(:,1),dA3deta))'*dRdeta - ...
    cross(dACovariant1(:,1),A3)'*ddRdetadeta;

%% Compute the variations of the scalar W_{alpha beta}, see equation (5e) - Nitsche Computation Vol. 10

% W11,U
W11KommaU = cross(dACovariant1(:,2),dACovariant2(:,1))'*dRdxi - ...
    cross(dACovariant1(:,2),dACovariant1(:,1))'*dRdeta;

% W22,U
W22KommaU = cross(dACovariant2(:,2),dACovariant2(:,1))'*dRdxi - ...
    cross(dACovariant2(:,2),dACovariant1(:,1))'*dRdeta;

% W12,U
W12KommaU = cross(dACovariant1(:,4),dACovariant2(:,1))'*dRdxi - ...
    cross(dACovariant1(:,4),dACovariant1(:,1))'*dRdeta;

%% Compute the derivatives of the variations of the scalar W_{alpha beta}, namely W_{alpha beta,gamma}, see equation (5f) - Nitsche Computation Vol. 10

% W11,1,U
W11Komma1KommaU = (cross(dACovariant1(:,3),dACovariant2(:,1)) + cross(dACovariant1(:,2),dACovariant1(:,4)))'*dRdxi + ...
    cross(dACovariant1(:,2),dACovariant2(:,1))'*ddRdxidxi - ...
    (cross(dACovariant1(:,3),dACovariant1(:,1)) + cross(dACovariant1(:,2),dACovariant1(:,2)))'*dRdeta - ...
    cross(dACovariant1(:,2),dACovariant1(:,1))'*ddRdxideta;

% W11,2,U
W11Komma2KommaU = (cross(dACovariant1(:,5),dACovariant2(:,1)) + cross(dACovariant1(:,2),dACovariant2(:,2)))'*dRdxi + ...
    cross(dACovariant1(:,2),dACovariant2(:,1))'*ddRdxideta - ...
    (cross(dACovariant1(:,5),dACovariant1(:,1)) + cross(dACovariant1(:,2),dACovariant1(:,4)))'*dRdeta - ...
    cross(dACovariant1(:,2),dACovariant1(:,1))'*ddRdetadeta;

% W22,1,U
W22Komma1KommaU = (cross(dACovariant1(:,6),dACovariant2(:,2)) + cross(dACovariant2(:,2),dACovariant1(:,4)))'*dRdxi + ...
    cross(dACovariant2(:,2),dACovariant2(:,1))'*ddRdxidxi - ...
    (cross(dACovariant1(:,6),dACovariant1(:,1)) + cross(dACovariant2(:,2),dACovariant1(:,2)))'*dRdeta - ...
    cross(dACovariant2(:,2),dACovariant1(:,1))'*ddRdxideta;

% W22,2,U
W22Komma2KommaU = (cross(dACovariant2(:,3),dACovariant2(:,2)) + cross(dACovariant2(:,2),dACovariant2(:,2)))'*dRdxi + ...
    cross(dACovariant2(:,2),dACovariant2(:,1))'*ddRdxideta - ...
    (cross(dACovariant2(:,3),dACovariant1(:,1)) + cross(dACovariant2(:,2),dACovariant1(:,4)))'*dRdeta - ...
    cross(dACovariant2(:,2),dACovariant1(:,1))'*ddRdetadeta;

% W12,1,U
W12Komma1KommaU = (cross(dACovariant1(:,5),dACovariant2(:,1)) + cross(dACovariant1(:,4),dACovariant1(:,4)))'*dRdxi + ...
    cross(dACovariant1(:,4),dACovariant2(:,1))'*ddRdxidxi - ...
    (cross(dACovariant1(:,5),dACovariant1(:,1)) + cross(dACovariant1(:,4),dACovariant1(:,2)))'*dRdeta - ...
    cross(dACovariant1(:,4),dACovariant1(:,1))'*ddRdxideta;

% W12,2,U
W12Komma2KommaU = (cross(dACovariant1(:,6),dACovariant2(:,1)) + cross(dACovariant1(:,4),dACovariant2(:,2)))'*dRdxi + ...
    cross(dACovariant1(:,4),dACovariant2(:,1))'*ddRdxideta - ...
    (cross(dACovariant1(:,6),dACovariant1(:,1)) + cross(dACovariant1(:,4),dACovariant1(:,4)))'*dRdeta - ...
    cross(dACovariant1(:,4),dACovariant1(:,1))'*ddRdetadeta;

%% 2. Compute the D matrix (see Nitsche - Computation Vol. 4) which is the variation of W vector (see equation (1c) Nitsche - Computation Vol. 5) 
D = - A2Cross*dRdxi + A1Cross*dRdeta;

%% 3. Compute the variation of the C vector (see Nitsche - Computation Vol. 4)
CVariation = 1/norm(A3Tilde)*D - 1/norm(A3Tilde)^3*A3Tilde*(A3Tilde'*D);

%% 4. Compute variations of vectors W,gamma (see equation (8) Nitsche - Computation Vol. 5)

% dWdxiVariation (see equation (20) Nitsche - Computation Vol. 5)
dWdxiVariation = - dA2dxiCross*dRdxi - A2Cross*ddRdxidxi + ...
    dA1dxiCross*dRdeta + A1Cross*ddRdxideta;

% dWdetaVariation (see equation (22) Nitsche - Computation Vol. 5)
dWdetaVariation = - dA2detaCross*dRdxi - A2Cross*ddRdxideta + ...
    dA2dxiCross*dRdeta + A1Cross*ddRdetadeta;

%% 5. Compute the variations of the derivatives of the strain components in the contravariant basis (see equation (2) Nitsche - Computation Vol. 5)

% dkappa11dxiVariation (see equation (6b) Nitsche - Computation Vol.
% 5 (Repetition))
dkappa11dxiVariation = - dA3dxi'*ddRdxidxi - A3'*dddRdxidxidxi + ...
    1/norm(A3Tilde)*(1/norm(A3Tilde)*dAKommaAlpha(1,1)*(dACovariant1(:,2) + ...
    BV(1,1)*A3) - (dACovariant1(:,3) + dBV11dxi*A3 + BV(1,1)*dA3dxi))'*...
    (-A2Cross*dRdxi + A1Cross*dRdeta) - ...
    1/norm(A3Tilde)*(dACovariant1(:,2) + BV(1,1)*A3)'*dWdxiVariation;

% dkappa11dxiKommaU, see equation (6) Nitsche - Computation Vol. 
% 10
dkappa11dxiKommaU = - A3'*dddRdxidxidxi - dA3dxi'*ddRdxidxi + ...
    (djInvdxi*BV(1,1) + jInv*dBV11dxi)*CKommaU + jInv*BV(1,1)*CKomma1KommaU + ...
    djInvdxi*W11KommaU + jInv*W11Komma1KommaU;

% dkappa11dxiVariation (see equation (6c) Nitsche - Computation Vol.
% 5 (Repetition))
dkappa11detaVariation = - dA3deta'*ddRdxidxi - A3'*dddRdxidxideta + ...
    1/norm(A3Tilde)*(1/norm(A3Tilde)*dAKommaAlpha(2,1)*(dACovariant1(:,2) + ...
    BV(1,1)*A3) - (dACovariant1(:,5) + dBV11deta*A3 + BV(1,1)*dA3deta))'*...
    (-A2Cross*dRdxi + A1Cross*dRdeta) - ...
    1/norm(A3Tilde)*(dACovariant1(:,2) + BV(1,1)*A3)'*dWdetaVariation;

% dkappa11detaKommaU, see equation (6) Nitsche - Computation Vol. 
% 10
dkappa11detaKommaU = - A3'*dddRdxidxideta - dA3deta'*ddRdxidxi + ...
    (djInvdeta*BV(1,1) + jInv*dBV11deta)*CKommaU + jInv*BV(1,1)*CKomma2KommaU + ...
    djInvdeta*W11KommaU + jInv*W11Komma2KommaU;

% dkappa22dxiVariation (see equation (11a) Nitsche - Computation Vol.
% 5 (Repetition))
dkappa22dxiVariation = - dA3dxi'*ddRdetadeta - A3'*dddRdetadetadxi + ...
    1/norm(A3Tilde)*(1/norm(A3Tilde)*dAKommaAlpha(1,1)*(dACovariant2(:,2) + ...
    BV(2,1)*A3) - (dACovariant1(:,6) + dBV22dxi*A3 + BV(2,1)*dA3dxi))'*...
    (- A2Cross*dRdxi + A1Cross*dRdeta) - ...
    1/norm(A3Tilde)*(dACovariant2(:,2) + BV(2,1)*A3)'*dWdxiVariation;

% dkappa22dxiKommaU, see equation (6) Nitsche - Computation Vol. 
% 10
dkappa22dxiKommaU = - A3'*dddRdetadetadxi - dA3dxi'*ddRdetadeta + ...
    (djInvdxi*BV(2,1) + jInv*dBV22dxi)*CKommaU + jInv*BV(2,1)*CKomma1KommaU + ...
    djInvdxi*W22KommaU + jInv*W22Komma1KommaU;

% dkappa11dxiVariation (see equation (11b) Nitsche - Computation Vol.
% 5 (Repetition))
dkappa22detaVariation = - dA3deta'*ddRdetadeta - A3'*dddRdetadetadeta + ...
    1/norm(A3Tilde)*(1/norm(A3Tilde)*dAKommaAlpha(2,1)*(dACovariant2(:,2) + ...
    BV(2,1)*A3) - (dACovariant2(:,3) + dBV22deta*A3 + BV(2,1)*dA3deta))'*...
    (- A2Cross*dRdxi + A1Cross*dRdeta) - ...
    1/norm(A3Tilde)*(dACovariant2(:,2) + BV(2,1)*A3)'*dWdetaVariation;

% dkappa22detaKommaU, see equation (6) Nitsche - Computation Vol. 
% 10
dkappa22detaKommaU = - A3'*dddRdetadetadeta - dA3deta'*ddRdetadeta + ...
    (djInvdeta*BV(2,1) + jInv*dBV22deta)*CKommaU + jInv*BV(2,1)*CKomma2KommaU + ...
    djInvdeta*W22KommaU + jInv*W22Komma2KommaU;

% dkappa11dxiVariation (see equation (14a) Nitsche - Computation Vol.
% 5 (Repetition))
dkappa12dxiVariation = - dA3dxi'*ddRdxideta - A3'*dddRdxidxideta + ...
    1/norm(A3Tilde)*(1/norm(A3Tilde)*dAKommaAlpha(1,1)*(dACovariant1(:,4) + ...
    BV(3,1)*A3) - (dACovariant1(:,5) + dBV12dxi*A3 + BV(3,1)*dA3dxi))'*...
    (- A2Cross*dRdxi + A1Cross*dRdeta) - ...
    1/norm(A3Tilde)*(dACovariant1(:,4) + BV(3,1)*A3)'*dWdxiVariation;

% dkappa12dxiKommaU, see equation (6) Nitsche - Computation Vol. 
% 10
dkappa12dxiKommaU = - A3'*dddRdxidxideta - dA3dxi'*ddRdxideta + ...
    (djInvdxi*BV(3,1) + jInv*dBV12dxi)*CKommaU + jInv*BV(3,1)*CKomma1KommaU + ...
    djInvdxi*W12KommaU + jInv*W12Komma1KommaU;

% dkappa11dxiVariation (see equation (14b) Nitsche - Computation Vol.
% 5 (Repetition))
dkappa12detaVariation = - dA3deta'*ddRdxideta - A3'*dddRdetadetadxi + ...
    1/norm(A3Tilde)*(1/norm(A3Tilde)*dAKommaAlpha(2,1)*(dACovariant1(:,4) + ...
    BV(3,1)*A3) - (dACovariant1(:,6) + dBV12deta*A3 + BV(3,1)*dA3deta))'*...
    (- A2Cross*dRdxi + A1Cross*dRdeta) - ...
    1/norm(A3Tilde)*(dACovariant1(:,4) + BV(3,1)*A3)'*dWdetaVariation;

% dkappa12dxiKommaU, see equation (6) Nitsche - Computation Vol. 
% 10
dkappa12detaKommaU = - A3'*dddRdetadetadxi - dA3deta'*ddRdxideta + ...
    (djInvdeta*BV(3,1) + jInv*dBV12deta)*CKommaU + jInv*BV(3,1)*CKomma2KommaU + ...
    djInvdeta*W12KommaU + jInv*W12Komma2KommaU;

%% 6. Compute the B-operator matrices for the membrane strain, the bending strain as well as the parametric derivatives of the bending strain in the covariant basis

% Compute the complete B-operator matrix for the membrane strain (see 
% Nitsche - Computation Vol. 4 - Checkpoint)
BOperatorMembraneStrainContra = [dACovariant1(:,1)'*dRdxi
                                 dACovariant2(:,1)'*dRdeta
                                 .5*(dACovariant2(:,1)'*dRdxi + dACovariant1(:,1)'*dRdeta)];

% Compute the complete B-operator matrix for the bending strain (see
% Nitsche - Computation Vol. 4 - Checkpoint)
% strain11 -> dG1/dxi
% strain22 -> dG2/deta
% strain12 -> dG1/deta = dG2/dxi
BOperatorBendingStrainContra = [- (A3Tilde'/norm(A3Tilde)*ddRdxidxi + dACovariant1(:,2)'*CVariation)
                                - (A3Tilde'/norm(A3Tilde)*ddRdetadeta + dACovariant2(:,2)'*CVariation)
                                - (A3Tilde'/norm(A3Tilde)*ddRdxideta + dACovariant1(:,4)'*CVariation)];
                         
% Compute the B-operator matrix for the parametric derivatives of the
% bending strain (see Nitsche - Computation Vol. 5)

% Parametric derivative of the bending strain in the xi-direction
dOperatorBendingStraindxiContra = [dkappa11dxiKommaU
                                   dkappa22dxiKommaU
                                   dkappa12dxiKommaU];
                            
% Parametric derivative of the bending strain in the eta-direction
dOperatorBendingStraindetaContra = [dkappa11detaKommaU
                                    dkappa22detaKommaU
                                    dkappa12detaKommaU];

end