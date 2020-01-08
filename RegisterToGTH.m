% reference reconstruction to ground truth
function [Q,alpha,signo]=RegisterToGTH(Q,Qg)
Qx=Q(1,:);
Qy=Q(2,:);
Qz=Q(3,:);
px=Qg(1,:);
py=Qg(2,:);
pz=Qg(3,:);
signo=1;
alpha=(inv(Qx(:)'*Qx(:)+Qy(:)'*Qy(:)+Qz(:)'*Qz(:))*(Qx(:)'*px(:)+Qy(:)'*py(:)+Qz(:)'*pz(:)));
Qx=alpha.*Qx;Qy=alpha.*Qy;Qz=alpha.*Qz;
error1=sqrt((norm(pz-Qz)^2+norm(px-Qx)^2+norm(py-Qy)^2)/length(px));
error2=sqrt((norm(pz+Qz)^2+norm(px+Qx)^2+norm(py+Qy)^2)/length(px));
if(error2<error1)
    signo=-1;
Qx=-Qx;Qy=-Qy;Qz=-Qz;
end
% Qx = Qx-mean(Qx);
% Qy = Qy-mean(Qy);


Q(1,:)=Qx;
Q(2,:)=Qy;
Q(3,:)=Qz;
