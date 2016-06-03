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
function sinusoidal_y = load_sinusoidal_y(p,i,u,U,q,j,v,V,CP)
%% Function documentation
%
% Returns Neumann boundary condition that has been recovered from 
% analytical solution for the recangular plate fixed at one edge and
% subject to a sinusoidal loading at the opposite edge
%
%  Input : 
%    p,q : Polynomial degrees
%    i,j : Knot span indeces
%    u,v : Parametric coordinates on the surface
%     CP : The set of control points and weights
%
% Output :
%  shear_y : The value of the analytical load at the surface location (u,v) 

%% Function main body

% The load maximum amplitude
Q = 1e3;

% The height of the plate
h = 2;

% Initialize the coordinates
% x = 0;
y = 0;

% Compute the basis functions affecting the knot span
Rb = nurbs_basis_functions2D(i,p,u,U,j,q,v,V,CP);

% initialize counter
k = 0;

for c = 0:q 
    for b = 0:p
        % update counter
        k = k + 1;
        
        % compute the location on x-y plane
        y = Rb(k)*CP(i-p+b,j-q+c,2) + y;  
    end
end

% Compute the load at the point x-y
% sinusoidal_y = Q*y/h - (3e-1)*Q*sin(pi*y/h); %Q*sin(pi*y/h)*cos(pi*y/h);
pot=1;
sinusoidal_y = 1e3*Q*(y^pot)/(h^pot);%cos(pi*y/h);

end

