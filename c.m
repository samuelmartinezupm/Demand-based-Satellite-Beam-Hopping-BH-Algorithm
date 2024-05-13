classdef c < handle
   properties (Access=public)
      number
      center
      radius
      users=[]
      % Computed
      betta_to_sat
      interfering
      % Dynamic (t=1:frame)
      aggregated_traffic=[]
      BW=[0]
      P=[0]
      active=[0]
      colour=[0]
   end
   properties (Constant)
   colour_colours=[{'blue'}, {'red'}, {'green'}, {'magenta'},{'cyan'}, {'black'}, {'yellow'}]
   end
   methods 
       function obj = c(number,center,radius,interfering) %constructor 
         obj.number = number;
         obj.center = center;
         obj.radius = radius;   
         obj.interfering = interfering;   
       end
       function pgon=draw(obj,t)
           pgon = nsidedpoly(1000, 'Center',  obj.center, 'Radius', obj.radius);
           [xc,yc] = centroid(pgon);
           if obj.colour(t)==0
            plot(pgon, 'FaceColor', char("white"))
           else
            plot(pgon, 'FaceColor', char(obj.colour_colours(obj.colour(t))))   
           end
           axis equal
           text(xc,yc, sprintf(num2str(obj.number), xc, yc), 'HorizontalAlignment','center', 'VerticalAlignment','middle','FontSize',14)
       end
       function pgon=draw_latlong(obj,lat_sp,lon_sp,t)
           pgon = nsidedpoly(1000, 'Center',  obj.center, 'Radius', obj.radius);
           Latitude=lat_sp+(pgon.Vertices(:,2))/110.574; %approximated conversion from km diff. in x and y to degrees in lat. and long.
           Longitude=lon_sp+(pgon.Vertices(:,1))./(111.32.*cosd(Latitude));
           plot(polyshape(Longitude',Latitude'),'FaceColor','none','EdgeColor', [mod(obj.active(t)+1,2),obj.active(t),0])
       end
       function obj=adduser(obj,u)
           obj.users=[obj.users u];
       end
       function obj=aggregatetraffic(obj)  % FALTA RESTA DE LO YA SERVIDO!
           sum=0;
           for u=1:length(obj.users)
               sum=sum+obj.users(u).traffic_pending;
           end
           obj.aggregated_traffic=[obj.aggregated_traffic sum];
       end
       function obj = compute_betta_to_sat(obj, h)  % Betta to Satellite from cell center (0,0):
           Re=6378; % Earth Radius [km]
           %1)
           gamma=(sqrt(obj.center(1)^2+obj.center(2)^2)*360)/(2*pi*Re);
           %2)
           h_prime=Re*(1-cos(gamma*pi/180));
           z=sin(gamma*pi/180)*Re;
           %3)
           betta=(atan(z/(h_prime+h)))*180/pi;
           obj.betta_to_sat=betta;
       end
   end
end

