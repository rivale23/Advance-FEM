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
function plot_ControlPolygonBSplineCurve(CP)
%% Function documenation
% 
% Plots Control Points and Control Polygon corresponding to a B-Spline 
% curve
%
%   Input :
%      CP : The control points for a B-Spline curve
%
%  Output :
%           graphics
%
%% Function main body

% Plot the intermediate Control Points and their connections
for l=1:length(CP(:,1))-1
    plot(CP(l,1),CP(l,2),'--or');
    plot(CP(l:l+1,1),CP(l:l+1,2),'--or');
end

end