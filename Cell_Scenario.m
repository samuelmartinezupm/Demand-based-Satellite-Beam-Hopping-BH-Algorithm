function [l,h,theta,D_footprint,number_cells,centers, interfering]=Cell_Scenario(h_sat,el_min,rings)

global display;

%% Data:
Re=6371; % Earth Radius [km]

%% Trigonometric raltionships:
betta_max=(asin((Re/(Re+h_sat))*cos(el_min*pi/180)))*180/pi;
gamma_max=90-betta_max-el_min;
h_prime_max=Re*(1-cos(gamma_max*pi/180));
A_footprint=2*pi*Re*h_prime_max;
D_footprint=((2*gamma_max)/360)*(2*pi*Re);

%% Cell Scenario - Dimensions [l,h]:
if mod(rings,2) % ODD RINGS -> radious sequence: an=9n^2-15n+7
    l=(D_footprint/2)/sqrt(9*((rings+1)/2)^2-15*((rings+1)/2)+7);
else % EVEN RINGS
    l=(D_footprint/2)/(rings+(rings/2-1));
end
h=(sqrt(3)/2)*l;
%1)
gamma=(2*l*360)/(2*pi*Re*2);
%2)
h_prime=Re*(1-cos(gamma*pi/180));
z=sin(gamma*pi/180)*Re;
%3)
betta=(atan(z/(h_prime+h_sat)))*180/pi;
%4)
theta=2*betta %FOV

% dmax -> Account for losses between [h<->dmax]
d_max=(Re+h_sat)*sqrt((Re/(Re+h_sat))^2+1-2*(Re/(Re+h_sat))*cos(gamma_max*pi/180));

%% Cell Scenario - Centers 
c_scenario=[];
number_cells=0;

for i=1:rings
    number_cells=number_cells+6*(i-1);
end
number_cells=number_cells+1;

centers=zeros(number_cells, 2);

cont=2; % Central one first index

for r=2:rings % Fill centers matrix by going into ring rounds, with 6 loops, each for each side of the hexagon.
    centers(cont,:)=[0,2*(r-1)*h]; %init (upper cell)
    cont=cont+1;
    for i=1:r-1
        centers(cont,:)=[centers(cont-1,1)+(l+l/2),centers(cont-1,2)-(h)];
        cont=cont+1;
    end
    for i=1:r-1
        centers(cont,:)=[centers(cont-1,1),centers(cont-1,2)-(2*h)];
        cont=cont+1;
    end
    for i=1:r-1
        centers(cont,:)=[centers(cont-1,1)-(l+l/2),centers(cont-1,2)-(h)];
        cont=cont+1;
    end
    for i=1:r-1
        centers(cont,:)=[centers(cont-1,1)-(l+l/2),centers(cont-1,2)+(h)];
        cont=cont+1;
    end
    for i=1:r-1
        centers(cont,:)=[centers(cont-1,1),centers(cont-1,2)+(2*h)];
        cont=cont+1;
    end
    for i=1:r-2
        centers(cont,:)=[centers(cont-1,1)+(l+l/2),centers(cont-1,2)+(h)];
        cont=cont+1;
    end
end

%% Interfering Cells Computation
d_interfering=2*h;
interfering=[];
for i=1:number_cells
    cont=0;
    for j=1:number_cells
        if i~=j
            if abs(sqrt((centers(i,1)-centers(j,1))^2+(centers(i,2)-centers(j,2))^2)-d_interfering)<1
                cont=cont+1;
                interfering(i,cont)=j;
            end
        end
    end
end

    %     % VISUALIZATION
    %     % Satellite Footprint:
    %     figure
    %     pgon = nsidedpoly(1000, 'Center',  [0,0], 'Radius', D_footprint/2);
    %     plot(pgon, 'FaceColor', 'y');
    %     axis equal
    %     hold on
    % 
    %     % Cell Footprint 
    %     for i=1:number_cells
    %         c_scenario=[c_scenario c(i,centers(i,:),l,1)];
    %         c_scenario(i).draw(1);
    %         hold on
    %     end
    %     title('Cell Scenario')
    %     hold off

end