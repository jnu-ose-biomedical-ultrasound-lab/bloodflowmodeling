

%--------------------------------------------------------------------------
 % rbcmodeling.m

 % Last updated: March 2019, LEE Cheong-Ah
 
 % Jeju National University-Biomedical Ultrasound Lab
 
 % Details: This is the primary function to run simulations. 
 %          Simulate time-dependent particle motion, and RBC aggregation under sinusoidal pulsatile flow. 

 % If you use our code, please cite our paper:
 % LEE, Cheong-Ah; KONG, Qi; PAENG, Dong-Guk. Depletion-model-based numerical simulation of the kinetics of red blood cell aggregation under sinusoidal pulsatile flow. Biorheology, 2018, Preprint: 1-13.
 
%--------------------------------------------------------------------------

%% particle motion dependent on time and RBC aggregation under sinusoidal pulsatile flow 
% this simulation code calculate the particle aggregation in the ROI
% under the sinusoidal pulsatile flow
% 2019.03.17
clc;clear all;close all;

currentFolder = pwd

%cd 'C:\Users\USL\Documents\MATLAB\RBC modeling code'
%% 1. Boundary of the tube
Bound_y=0.1e-3;                                      % Tube diameter [m]
Bound_x=0.3e-3;                                      % Tude length [m]
P_R_x(1,:)  = linspace(0,Bound_x,250);               % grid of x
P_R_y(:,1)  = linspace(-Bound_y/2,Bound_y/2,250);    % grid of y
full_area = Bound_x*Bound_y; 
%% 2. Set the input parameters
mass=9.8*10^(-14); % Mass of RBC [kg]
R=4*10^(-6);       % Radius of RBC [m]
Vs_R=2.7e-6;       % blood viscosity [Pa.s]
PNUM= floor(n_RBC(Bound_x,Bound_y,R,0.4)); % hematocrit (the number of particle)
AGG=2;             % Choice the RBC aggregation parameter
%% 3. Flow condition
str_rate=60;                % Stroke rate [bpm]
freq=str_rate/60;           % Fluid frequency with stroke rate [Hz]
Vp_flow=0.3e-3;             % Wave length [m]                                              
wave_angle=freq*(2*pi);     % 2*pi*f [rad/sec]
wave_number=2*pi/Vp_flow;   % Wavenumber=2*pi/wave length [rad/m]
in_t=0;                     % Phase[radian] 
AC_R=zeros(PNUM,2);         % Intial particle acceleration vector
dv=0.0015;                  % Velocity amplitude [m/s]
v_mean=0.0035;              % Mean velocity [m/s]
period =2;                  % period
%% 4. RBC raw position
V_R=zeros([PNUM 2]);        %particle velocity field                                      
lo=dist_RBC(Bound_x-2*R,Bound_y-2*R,PNUM,R);    % Distance of the particles (x,y)
P_R=[lo(:,1)+R lo(:,2)];                        % Position of the Particles (x,y)
%% 5. make the movie
% aviobj=VideoWriter('RBC_motion.avi');
% aviobj.FrameRate = 1;                                                        
% open(aviobj);
 
%% 6. Main presses part
once_t=0; indx = 0; kkk=0;
for T=1:10000*period;
    t=[0.1e-3]; T_num = T*t;
    cd(currentFolder)
    %cd 'C:\Users\USL\Documents\MATLAB\RBC modeling code';
   %% boundary condition for location of particle x,y  
    Ovx_index1=find(P_R(:,1)>Bound_x);
    P_R(Ovx_index1,1)=P_R(Ovx_index1,1)-(Bound_x);
    Ovx_index2=find(P_R(:,1)<=0);
    P_R(Ovx_index2,1)=P_R(Ovx_index2,1)+(Bound_x);
    Ovx_indey1=find(P_R(:,2)>(Bound_y/2-R/2));
    P_R(Ovx_indey1,2)=P_R(Ovx_indey1,2)-Bound_y+R;
    Ovx_indey2=find(P_R(:,2)<-Bound_y/2+R/2);
    P_R(Ovx_indey2,2)=P_R(Ovx_indey2,2)+Bound_y-R;   
   %% This part calculate new position of the RBCs
    velocity_particle   = (dv*sin(wave_angle*(T_num)-wave_number*P_R(:,1)+in_t)+v_mean).*(1-(4.*(P_R(:,2).^2)./(Bound_y^2)));    
    velocity_field      = (dv*sin(wave_angle*(T_num)-wave_number*P_R_x(1,:)+in_t)+v_mean)'*(1-(4.*(P_R_y(:,1).^2)./(Bound_y^2)))';   % velocity field [m/s]
   %% calculate net force(hydro dynamic force, aggregation force, elastic force)
    uu= velocity_particle';                                    % particle's velocity
    [d_R FA_un d_RF]=GetMD(P_R,Bound_x,R);                     % Get the distance and direction among particles
    FA_R = GetMF(d_R,AGG,R).*FA_un;    
    Fh_R=6*pi*Vs_R*R*([uu',zeros(PNUM,1)]-V_R);                % Hydrodynamic force
    F_Rt=[real(sum(FA_R, 1))' imag(sum(FA_R, 1))'];            % Calculation of the interactional forces
    F_R=[real(sum(FA_R, 1))' imag(sum(FA_R, 1))'] +Fh_R;       % Total force of particle
    AC_R=F_R./mass;                                            % Acceleration of the particles(Newton's second law)
    Step_xy=V_R.*t+0.5*AC_R.*t^2;                              % Positon(x,y) of the particles
    V_R=V_R+AC_R.*t;                                         % velocity function [m/s]                     
    P_R=P_R+V_R.*t+0.5*AC_R.*t^2;                            % particle position [m]  
    once_t=once_t+t;                                         % time step [sec]   
    kkk=kkk+1;      
   %% the part indicate the velocity, acceleration, shear rate in the ROI rage 
    acceleration =  (dv*wave_angle*cos(wave_angle*(T_num)-wave_number*P_R_x(1,:)+in_t))'*(1-(4.*(P_R_y(:,1).^2)./(Bound_y^2)))';  % acceleration field [m/s2]
    win_acc_data = unique(acceleration); win_acc_data = sort(win_acc_data); 
    Acc_mean = mean(win_acc_data);
    shearrate = ((dv*sin(wave_angle*(T_num)-wave_number*P_R_x(1,:)+in_t)+v_mean)'*((8.*abs(P_R_y(:,1))./((Bound_y)^2)))');        % shear rate in the y direction
%   shearrate_x = (-dv*wave_number*cos(wave_angle*(T_num)-wave_number*P_R_x(1,:)+in_t)).*(1-(4.*(P_R_y(:,1).^2)./(Bound_y^2)));  % shear rate in the x direction
%   shearrate_magnitude = sqrt(shearrate_x.^2+shearrate_y.^2);                                                                   % magnitude field of shear rate 
    win_shear_data = unique(shearrate); win_shear_data = sort(win_shear_data);
    SR_mean = mean(win_shear_data);   
 
    if  once_t >=5e-3
    %% Analysis of the aggregate
        [d_R FA_un d_RF]=GetMD(P_R,Bound_x,R); 
        AZ=d_RF<2.02*R;                                      
        AZ1=AZ;
        AZ=triu(AZ);                                                
        [AGG_x,AGG_y]=find(AZ==1);                   
        AGG_p=[AGG_x AGG_y]; 
        AGG_RBC_No=unique(AGG_p); 
        AGG_RBC_No2 = AGG_RBC_No; % aggregated RBC number
        [remove] = find(P_R(AGG_RBC_No,2)>(Bound_y/2-2*R)|P_R(AGG_RBC_No,2)<(-Bound_y/2+2*R));  
        AGG_RBC_No(remove)=[];    
    %% Calculation the velocity, acceleration, shear rate of the each particles
        particle_acc = (dv*wave_angle*cos(wave_angle*(T_num)-wave_number*P_R(AGG_RBC_No,1)+in_t)).*(1-(4.*(P_R(AGG_RBC_No,2).^2)./(Bound_y^2)));
        particle_sh = ((dv*sin(wave_angle*(T_num)-wave_number*P_R(AGG_RBC_No,1)+in_t)+v_mean).*((8.*abs(P_R(AGG_RBC_No,2))./((Bound_y)^2))));    
    %% area calculation system                     
        idcs = strfind(currentFolder,'\');
        newdir = currentFolder(1:idcs(end)-1);
        cd(newdir)
 %   cd currentFolder+'\test' %  dir:  new folder for the save the data
        save(['data_' num2str(indx)  '.mat'],'v_mean','dv','R', 'Bound_y','Bound_x','T_num','full_area','win_acc_data','SR_mean','acceleration','particle_acc','Acc_mean','AGG_RBC_No','AGG_RBC_No2','P_R','win_shear_data','shearrate','particle_sh','P_R_x','P_R_y')
        indx = indx+1;
%% check the figure %%    
    cd(currentFolder)
%cd 'C:\Users\USL\Documents\MATLAB\RBC modeling code'
    figure(1)
     imge1 = imagesc(P_R_x, P_R_y,shearrate)
     colormap(jet)  
     colorbar ();
     hold on 
     set(gcf,'position', [0 169 1200 515])
     DrawCircle(P_R(:,1),P_R(:,2),R, 100, 'r.')
     axis equal
     hold on
     DrawCircle(P_R(AGG_RBC_No,1),P_R(AGG_RBC_No,2),R, 100, 'k.')
     barh(-Bound_y/1.9,Bound_x+2*R,1e-6);
     barh(Bound_y/1.9,Bound_x+2*R,1e-6);
     hold off
     title(['Time: ',num2str(T*t),'s'])
     xlim([0 Bound_x])
     ylim([-Bound_y/2-2e-5 Bound_y/2+2e-5])
%      M=getframe(figure(1));
%      writeVideo(aviobj,M)
        once_t = 0; 
    end
 
end
%  
