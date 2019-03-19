
%--------------------------------------------------------------------------
 % GetMD.m

 % Last updated: March 2019, LEE Cheong-Ah
 
 % Jeju National University-Biomedical Ultrasound Lab
 
 % Details: Get the distance and direction among particles. 

 % If you use our code, please cite our paper:
 % LEE, Cheong-Ah; KONG, Qi; PAENG, Dong-Guk. Depletion-model-based numerical simulation of the kinetics of red blood cell aggregation under sinusoidal pulsatile flow. Biorheology, 2018, Preprint: 1-13.
 
%--------------------------------------------------------------------------



% Get the distance and direction among particles（计算RBC的距离和方向）
% Rnm is the matrix of distances among particles
%Unit_V is the unit direction vector among particles
%Arrange like
%Rnm（d_R）:  表示距离：第一列代表第一个点到其他点的距离(不包括本身的点)，第二列代表第二个点到其他点的距离
%Unit_V: 表示距离的复数部分，x+yi，x是横坐标，y是纵坐标，abs（Unit_V）=Rnm
%d_RF:  表示距离：第一列代表第一个点到其他点的距离(包括本身的点)，第二列代表第二个点到其他点的距离
              %中间的对角线本来代表点到点本身的距离（0），但是这里变成100,

function [Rnm Unit_V d_RF]=GetMD(P_R,x_boundy,R)
P_n=size(P_R,1);                                                  %多少个RBC数
One_Pn=ones([1,P_n]);                                        %都是单位1的一行
One_Pn_1=ones([1,P_n-1]);                                 %RBC总数减一

%xn_xm
P_x=P_R(:,1);                                                        %第一列，11个点的横坐标
P_xmn=(P_x*One_Pn);                                          %每一列都相同的方阵，横坐标都铺开？
P_xmn_1=(P_x*One_Pn_1)';                                  %每一行都相同的十行，（少了一列）
b=triu(P_xmn,1);                                                  %从第二行到最后一行
b=b(1:end-1,:);                                                     %一列算是从第一个到最后一个点的横坐标，从第一个点到所有点
c=tril(P_xmn,-1);
c=c(2:end,:);                                                         % 第一列是所以点的横坐标倒过来， 从所有点到第一个点
P_xx=b+c;                                                            %奇怪的组合
P_xn_xm=P_xmn_1-P_xx;                                      %第一列代表第一个点到其他点的距离，第二列代表第二个点到其他点的距离
                                                                             %以此类推(但只是横坐标)
%yn_ym
P_y=P_R(:,2);
P_ymn=(P_y*One_Pn);
P_ymn_1=(P_y*One_Pn_1)';
b=triu(P_ymn,1);
b=b(1:end-1,:);
c=tril(P_ymn,-1);
c=c(2:end,:);
P_yy=b+c;
P_yn_ym=P_ymn_1-P_yy;
P_xn_xm(P_xn_xm>x_boundy-4*R)=P_xn_xm(P_xn_xm>x_boundy-4*R)-x_boundy;
P_xn_xm(P_xn_xm<-(x_boundy-4*R))=P_xn_xm(P_xn_xm<-(x_boundy-4*R))+x_boundy;

%Rnm
Rnm=sqrt(P_xn_xm.^2+P_yn_ym.^2);                %表示距离：第一列代表第一个点到其他点的距离，第二列代表第二个点到其他点的距离
%Angle
Unit_V=(P_xn_xm+P_yn_ym*i)./Rnm;

%d_RF full size D to D
P_xf=P_xmn-P_xmn';                                           %表示横坐标：第一行代表第一个点到所有点的距离（包括本身的点）
P_xf(P_xf>x_boundy-4*R)=P_xf(P_xf>x_boundy-4*R)-x_boundy;
P_xf(P_xf<-(x_boundy-4*R))=P_xf(P_xf<-(x_boundy-4*R))+x_boundy;
    %第二行代表第二个点到所有点的距离（包括本身的点）
P_yf=P_ymn-P_ymn';
d_RF=sqrt(P_xf.^2+P_yf.^2);                               %表示总的距离
d_RF=d_RF+100*diag(One_Pn);                          %对角线为什么加100？？？
