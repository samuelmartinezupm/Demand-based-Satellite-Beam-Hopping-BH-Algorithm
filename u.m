classdef u < handle
   properties (Access=public)
      % Static
      location
      type_station
      traffic_demand
      gain
      T_noise
      % Computed
      betta_to_sat
      distace_to_sat
      elevation_to_sat
      betta_to_center
      % Dynamic (t=1:frame)
      C_N=[]
      P=[]
      BW=[]
      I
      traffic_served=[]
      traffic_pending=0
   end
   methods
       function obj = u(location,type_station,traffic_demand,gain,T_noise) %constructor
         obj.location = location;
         obj.type_station = type_station;
         obj.traffic_demand = traffic_demand;
         obj.gain=gain;
         obj.T_noise=T_noise;
       end
       function obj = compute_distance_elevation_betta_to_sat(obj, h)  % Distance and Elevation to Satellite (0,0):
           Re=6378; % Earth Radius [km]
           %1)
           gamma=(sqrt(obj.location(1)^2+obj.location(2)^2)*360)/(2*pi*Re);
           %2)
           h_prime=Re*(1-cos(gamma*pi/180));
           z=sin(gamma*pi/180)*Re;
           %3)
           betta=(atan(z/(h_prime+h)))*180/pi;
           obj.betta_to_sat=betta;
           %4)
           obj.distace_to_sat=(h+h_prime)/cos(betta*pi/180);
           obj.elevation_to_sat=90-betta-gamma;
       end
       function obj = compute_betta_to_cell_center(obj, h, c)  % Betta computation as if center of the cell would have been the user, then compare the difference with respect to the real user and evaluate the real apperture from cell center.
           Re=6378; % Earth Radius [km]
           %1)
           gamma=(sqrt((obj.location(1)-c.center(1))^2+(obj.location(2)-c.center(2))^2)*360)/(2*pi*Re);
           %2)
           h_prime=Re*(1-cos(gamma*pi/180));
           z=sin(gamma*pi/180)*Re;
           %3)
           betta=(atan(z/(h_prime+h)))*180/pi;
           %4)
           obj.betta_to_center=betta; %abs(obj.betta_to_sat-betta)
       end
       function draw(obj)
       plot(obj.location(1),obj.location(2),'-s','MarkerSize',7,'MarkerEdgeColor','black','MarkerFaceColor','black') 
       end
   end
end

