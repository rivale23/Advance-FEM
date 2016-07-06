function [ SensitivityMatrix,Mass,Disp,BSplinePatch  ] = GetAllSensitivities(BSplinePatch, vectors, IndependentDirectionsFlag)
%GETALLSENSITIVITIES Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
   IndependentDirectionsFlag=false; 
end
[ SensitivityMatrix,Mass,Disp,BSplinePatch ] =Sensitivity_wrapper( BSplinePatch, vectors,'all',IndependentDirectionsFlag);


end
