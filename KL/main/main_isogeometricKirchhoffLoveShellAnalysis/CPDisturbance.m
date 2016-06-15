function [ BSplinePatchDist ] = CPDisturbance( BSplinePatch, Position, Direction, Magnitude,s )
%CPDISTURBANCE 
%gets the length of the vector of positions
count=size(Position,1);
X=zeros(1,count);
Y=X;
Z=X;
%gets the values of the components of the desired points
for i=1:count 
    xx=Position(i,1);
    yy=Position(i,2);
    X(i)=BSplinePatch.CP(xx,yy,1);  %gets X component
    Y(i)=BSplinePatch.CP(xx,yy,2); %gets y component
    Z(i)=BSplinePatch.CP(xx,yy,3);    %gets z component
end 

%gets the unit vector
Disturbance=Direction/norm(Direction);

%sets the magnitude
Disturbance=Disturbance.*Magnitude;

%creates a vector with all the modified CPs 
CPmod=zeros(3,count);
for i=1:count 
    CPmod(i,1)=X(i) +Disturbance(1);
    CPmod(i,2)=Y(i) +Disturbance(2);
    CPmod(i,3)=Z(i) +Disturbance(3);
end
CPs=BSplinePatch.CP;
for i=1:count 
    CPs(Position(i,1),Position(i,2),1)=CPmod(i,1);
    CPs(Position(i,1),Position(i,2),2)=CPmod(i,2);
    CPs(Position(i,1),Position(i,2),3)=CPmod(i,3);
end
BSplinePatchDist=BSplinePatch;

if(s==1)    
    BSplinePatchDist.CPd=CPs;
else
    BSplinePatchDist.CP=CPs;
end





end

