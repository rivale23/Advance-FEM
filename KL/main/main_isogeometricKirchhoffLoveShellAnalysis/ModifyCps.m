function [ ModifiedBspline ] = ModifyCps( BSplinepatch,DesiredCPs ,delta)

%   read all the desired CPs from the DesiredCPs variable and modify all of
%   them

%for test Ill just modify some
CP=BSplinepatch.CP;
CP(:,2,3)=CP(:,2,3)+delta;

ModifiedBspline=BSplinepatch;

ModifiedBspline.CP=CP;
end
