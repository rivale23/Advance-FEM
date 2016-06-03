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
function plot_BSplineCurveOnBSplineSurface...
    (p,q,Xi,Eta,CP,isNURBS,grid,xi1,eta1,xi2,eta2,color,lineWidth)
%% Function documentation
%
% Plots the parametric curve on the NURBS surface
%
%     Input :
%       p,q : Polynomial degrees
%    Xi,Eta : Knot vectors in u,v-direction
%        CP : Set of control points and weights
%      grid : Number of sampling points to be used
%  xi1,eta1 : Starting point parameters of the curve
%  xi2,eta2 : Ending point parameters of the curve
%     color : The color of the line to be plotted
% lineWidth : The width of the line
%
%   Output : graphics
%
% Function layout :
%
% 1. Get the coordinates of the sampling points on the curve
%
% 2. Create the geometry
%           
%% Function main body

%% 1. Get the coordinates of the sampling points on the curve
[Xp,Yp,Zp] = createBSplineCurveOnBSplineSurface...
    (p,q,Xi,Eta,CP,isNURBS,grid,xi1,eta1,xi2,eta2);

%% 2. Create the geometry
line(Xp,Yp,Zp,'Linewidth',lineWidth,'color',color);

end