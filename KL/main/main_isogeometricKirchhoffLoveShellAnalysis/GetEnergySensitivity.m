function [ SensitivityMatrix,SensitivityMass,BSplinePatch ] = GetEnergySensitivity( BSplinePatch, vectors, IndependentDirectionsFlag)

if nargin < 3
   IndependentDirectionsFlag=false; 
end
[ SensitivityMatrix,SensitivityMass,~,BSplinePatch ] =Sensitivity_wrapper( BSplinePatch, vectors,'energy',IndependentDirectionsFlag);

end