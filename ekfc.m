function [A, B, C, D]=ekfc(y, u, t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extended Kalman Filter Based Online System Identification and Denoising             %
%                                                                                     %
% Modify observation matrix (H), covariance matrix of the system noise (Q)            %  
% and covariance matrix of observation noise (R) for your model to get better results %
%                                                                                     %
% Author: Reza Samarghandi                                                            %
%                                                                                     %
% Email: Rezasamargandi@yahoo.com                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%y: Observed Parameters
%u: Control Inputs


C0=eye(n); %C Matrix of State Space (using this mean we call observe all of the states)
D0=zeros(n,m); %D Matrix of State Space (using this mean we don't use feed forward)
H=[C0 D0 zeros(n,(n^2+n*m))]; %observation matrix (using this mean we call observe all of the states and no feed-forward)
Q=1000.*(H.'*H); %covariance matrix of the system noise
R=1000.*eye(n,n); %covariance matrix of observation noise wk


n=size(y,1);
m=size(u,1);


if t==0
    y0=[y; u];    
    I0=eye(n,n+m);
    It0=I0.';
    I0v = It0(:);
    theta=[y0; I0v]; %state variable vector 
    gam=1000; %large value
    Pk=gam.*eye(n+m+n*(n+m)); %covariance matrix
end

Atp=theta(n+m+1:end,:);
Atil=reshape(Atp,[n+m n]).';

xtil=theta(1:n+m,:);

X=Atil*xtil;

X=[X; u];

Xh=X.';

bbc= kron(eye(n),Xh);

thetak1=[X; Atp];

F=[Atil, bbc;
   zeros(m+n*(n+m),n), eye(m+n*(n+m))];%Jacobian matrix

Pk1=F*Pk*F.'+Q;

S=H*Pk1*H.'+R; %noise covariance matrix

Kk=Pk1*H.'*S^-1; %Kalman gain


deltheta=Kk*(y-H*thetak1); %modification vector for state variables 

theta=thetak1+deltheta;

Pk=(eye(n+m+n*(n+m))-Kk*H)*Pk1;

Atp=theta(n+m+1:end,:);
Ap=reshape(Atp,[n+m n]).'; %[A B]


A=Ap(:,1:n); %A Matrix of State Space 
B=Ap(:,n+1:end); %B Matrix of State Space
C=H(:,1:n); %C Matrix of State Space
D=H(:,n+1:n+m); %D Matrix of State Space


end
