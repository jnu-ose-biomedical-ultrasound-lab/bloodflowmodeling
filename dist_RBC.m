function lo=dist_RBC(x,y,n,R)

%--------------------------------------------------------------------------
 % dist_RBC.m

 % Last updated: March 2019, LEE Cheong-Ah
 
 % Jeju National University-Biomedical Ultrasound Lab
 
 % Details: Generate random red blood cell (RBC) distribution.

 % If you use our code, please cite our paper:
 % LEE, Cheong-Ah; KONG, Qi; PAENG, Dong-Guk. Depletion-model-based numerical simulation of the kinetics of red blood cell aggregation under sinusoidal pulsatile flow. Biorheology, 2018, Preprint: 1-13.
 
%--------------------------------------------------------------------------


%RBC randomly distribution
lo=rand(n,2);
% x=x-4*R;
lo(:,1)=lo(:,1)*x;
lo(:,2)=(lo(:,2)-0.5)*y;
% r2=(2.7e-6*2)^2*1.1;
r2=(R*2)^2*1.1;
for k=1:n
    d=(lo(k,1)-lo(:,1)).^2+(lo(k,2)-lo(:,2)).^2;
    d(k)=inf;
    while min(d)<r2
        lo(k,:)=[rand()*x (rand()-0.5)*y];
        d=(lo(k,1)-lo(:,1)).^2+(lo(k,2)-lo(:,2)).^2;
        d(k)=inf;
    end
    
end
