function [Case,sg] = CaseSg(p1,p2,p3,PC,PD,v,u,n)
%This code serves for SubP and SubP_alpha to do the root-finding procedure


Case = zeros(n,6);
sg = inf*ones(n,2);
VC = (v + 1).^2;
VD = (v - 1).^2;
U = u.^2;


%case 7
I7 = find(p1==0);
sg(I7,[1 2]) = [abs(u(I7)) inf*ones(size(I7))];
Case(I7,[1 4 6]) = [ones(size(I7)) U(I7) U(I7)];

%Case 1
I11 = (p1>0).*(p3>0).*(u>0);
I1= find(I11>0);
Case(I1,[1 3 5]) = [VD(I1) VD(I1) VD(I1)];

% Case 2
I22 = (p1>0).*(p3>0).*(u<0);
I2= find(I22>0);
sg(I2,:) = [PC(I2) PD(I2)];
Case(I2,[1 4 5]) = [VC(I2) U(I2) VD(I2)];



%case 3
I33 = (p1>0).*(p3<0).*(u>0) + (p1<0).*(p2>0).*(u>0) + (p1<0).*(p2==0).*(u>0);
I3= find(I33>0);
sg(I3,[1 2]) = [PD(I3) inf*ones(size(I3))];
Case(I3, [1 4 6]) = [VD(I3) U(I3) U(I3)];


%case 4
I44 = (p1>0).*(p3<0).*(u<=0) + (p1>0).*(p3==0).*(u<0) + (p1<0).*(p2>0).*(u<=0);
I4= find(I44>0);
sg(I4,[1 2]) = [PC(I4) inf*ones(size(I4))];
Case(I4, [1 4 6]) = [VC(I4) U(I4) U(I4)];

% Case 5
I55 = (p1<0).*(p2<0).*(u>0);
I5= find(I55>0);
sg(I5,:) = [PD(I5) PC(I5)];
Case(I5,[1 4 5]) = [VD(I5) U(I5) VC(I5)];

% Case 6
I66 = (p1<0).*(p2<0).*(u<=0);
I6= find(I66>0);
Case(I6,[1 3 5]) = [VC(I6) VC(I6) VC(I6)];



end