function [ BSplinePatchDist ] = CPDisturbance( BSplinePatch, Position, Direction, Magnitude )
%CPDISTURBANCE Summary of this function goes here
%   Detailed explanation goes here
X=BSplinePatch.CP(Position(1),Position(2),1);  %gets X component
Y=BSplinePatch.CP(Position(1),Position(2),2); %gets y component
Z=BSplinePatch.CP(Position(1),Position(2),3);    %gets z component


Disturbance=Direction/norm(Direction);

Disturbance=Disturbance.*Magnitude;

CPmod=[X Y Z]+Disturbance;

CPs=BSplinePatch.CP;
CPs(Position(1),Position(2),1)=CPmod(1);
CPs(Position(1),Position(2),2)=CPmod(2);
CPs(Position(1),Position(2),3)=CPmod(3);

BSplinePatchDist.CPd=CPs;




end

