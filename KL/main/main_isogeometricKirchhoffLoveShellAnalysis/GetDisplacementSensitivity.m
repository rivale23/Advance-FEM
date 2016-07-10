function [ SensitivityDisplacement ,SensitivityMass,BSplinePatch] = GetDisplacementSensitivity( BSplinePatch, vectors, IndependentDirectionsFlag )

if nargin < 3
   IndependentDirectionsFlag=false; 
end
[ ~,SensitivityMass,SensitivityDisplacement ,BSplinePatch] =Sensitivity_wrapper( BSplinePatch, vectors,'displacement',IndependentDirectionsFlag);

end
