function [area Area_div win_para_data] = cjddk_pix_area(div, full_area,win_para_data,parameter,Bound_x, Bound_y )

%--------------------------------------------------------------------------
 % cjddk_pix_area.m

 % Last updated: March 2019, LEE Cheong-Ah
 
 % Jeju National University-Biomedical Ultrasound Lab
 
 % Details: Output the area based on divisions and input parameters. 

 % If you use our code, please cite our paper:
 % LEE, Cheong-Ah; KONG, Qi; PAENG, Dong-Guk. Depletion-model-based numerical simulation of the kinetics of red blood cell aggregation under sinusoidal pulsatile flow. Biorheology, 2018, Preprint: 1-13.
 
%--------------------------------------------------------------------------



%% input parameter
% div = 40;
% win_para_data  = win_acc_data;
% parameter = accelaration;
unit_area = full_area/div;
for i = 1:size(win_para_data)
        i
    parameter1 = [];
    parameter1 = parameter;
    parameter1(parameter1>=win_para_data(i) )=0;
%     figure(1)
%     imagesc(parameter1)
    
    pixel = numel(parameter1(parameter1~=0));
    dx = Bound_x/size(parameter1,2); dy = Bound_y/size(parameter1,1);      % unit = m
    area(:,i) = pixel*(dx*dy); % unit = m^2
end

    Area = 0; k = 1;
   [num]= find(area(1,:)>=unit_area);
   Area_div(:,k) = [area(1,num(1)); (full_area-area(1,num(1))); num(1)];
   
    for i = num(1):size(area,2)-1
        Area = Area + (area(:,i+1)-area(:,i));
        if Area>= unit_area
             k = k+1;
            Area_div(:,k) = [Area ;(full_area - Area);i];
            Area = 0;
        end
    end
    Area_div(:,k+1) =  [(full_area - area(1,Area_div(3,k))); 0; size(area,2)];
end
%   
   
