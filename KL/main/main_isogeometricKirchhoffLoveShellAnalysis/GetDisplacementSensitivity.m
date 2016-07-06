function [ SensitivityDisplacement ,b,BSplinePatch] = GetDisplacementSensitivity( BSplinePatch, vectors, IndependentDirectionsFlag )

if nargin < 3
   IndependentDirectionsFlag=false; 
end
[ a,b,SensitivityDisplacement ,BSplinePatch] =Sensitivity_wrapper( BSplinePatch, vectors,'displacement',IndependentDirectionsFlag);

end
