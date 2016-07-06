function [ SensitivityMatrix,a,BSplinePatch ] = GetEnergySensitivity( BSplinePatch, vectors, IndependentDirectionsFlag)

if nargin < 3
   IndependentDirectionsFlag=false; 
end
[ SensitivityMatrix,a,b,BSplinePatch ] =Sensitivity_wrapper( BSplinePatch, vectors,'energy',IndependentDirectionsFlag);

end