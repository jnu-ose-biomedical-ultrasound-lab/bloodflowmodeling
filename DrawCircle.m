function [BW]=DrawCircle(x, y, r, nseg, S)
% Draw a circle on the current figure using ploylines
%
%  DrawCircle(x, y, r, nseg, S,I)
%  A simple function for drawing a circle on graph.
%
%  INPUT: (x, y, r, nseg, S)
%  x, y:    Center of the circle
%  r:       Radius of the circle
%  nseg:    Number of segments for the circle
%  S:       Colors, plot symbols and line types
%
%  OUTPUT: None
%
%  BUG REPORT:
%  Please send your bug reports, comments and suggestions to
%  pengtao@glue.umd.edu . Thanks.

%  Author:  Tao Peng
%           Department of Mechanical Engineering
%           University of Maryland, College Park, Maryland 20742, USA
%           pengtao@glue.umd.edu
%  Version: alpha       Revision: Jan. 10, 2006


%--------------------------------------------------------------------------
 % DrawCircle.m

 % Last updated: March 2019, LEE Cheong-Ah
 
 % Jeju National University-Biomedical Ultrasound Lab
 
 % Details: Circle drawing function. Added to RBC blood flow sim. 

 % If you use our code, please cite our paper:
 % LEE, Cheong-Ah; KONG, Qi; PAENG, Dong-Guk. Depletion-model-based numerical simulation of the kinetics of red blood cell aggregation under sinusoidal pulsatile flow. Biorheology, 2018, Preprint: 1-13.
 
%--------------------------------------------------------------------------




size_x=size(x,1);
M_one=ones(1,size_x);
T_one=ones(1,nseg+1);
theta = 0 : (2 * pi / nseg) : (2 * pi);
pline_x = r *cos(theta)'*M_one + T_one'*x';
pline_y = r *sin(theta)'*M_one + T_one'*y';


fill(pline_x, pline_y, S);
%BW = roipoly(I,pline_x, pline_y);
