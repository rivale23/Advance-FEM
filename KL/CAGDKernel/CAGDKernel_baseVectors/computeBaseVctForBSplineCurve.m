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
function G = computeBaseVctForBSplineCurve(knotSpanIndex,p,CP,dR)
%% Function documentation
%
% Returns the base vector over a B-Spline curve.
%
%         Input :
% knotSpanIndex : The span where u is contained
%             p : The polynomial degree of the curve
%            CP : The Control points of the NURBS curve
%            dR : The basis functions and their derivatives at a given 
%                 parametric location
%
%        Output :
%             G : The base vector to the B-Spline curve at the given
%                 parametric location
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the base vector
%
%% Function main body

%% 0. Read input

% Tangent vector
G = zeros(3,1);

%% 1. Compute the base vector
for b = 0:p
    for a = 1:3
        % Compute the index
        index = knotSpanIndex - p + b;
        
        % Compute the tangent base vector to the curve
        G(a,1) = G(a,1) + dR(b+1,2)*CP(index,a);
    end
end

end