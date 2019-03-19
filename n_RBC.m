function n=n_RBC(x,y,r_RBC,con) 

%--------------------------------------------------------------------------
 % n_RBC.m

 % Last updated: March 2019, LEE Cheong-Ah
 
 % Jeju National University-Biomedical Ultrasound Lab
 
 % Details: Simulate time-dependent particle motion, and RBC aggregation under sinusoidal pulsatile flow. 

 % If you use our code, please cite our paper:
 % LEE, Cheong-Ah; KONG, Qi; PAENG, Dong-Guk. Depletion-model-based numerical simulation of the kinetics of red blood cell aggregation under sinusoidal pulsatile flow. Biorheology, 2018, Preprint: 1-13.
 
%--------------------------------------------------------------------------



%¼ÆËãRBC¸öÊý
s_RBC=x*y*con;
% r_RBC=2.7e-6;
s_RBC1=pi*r_RBC^2;
n=s_RBC/s_RBC1;
