clear all
clc
close all

global display;
display=0; % 0) No figures displayed: fast execution. 1) Display all figures.

% Fixed Parameters:
realization=1; %Realization Number (Monte Carlo)
colours=4; %Number of colours to consider. 1: Single Colour. 4: 2 freq. + 2 pol.
beams=8; % Number of simultaneous beams.
TTL=5; % Weighting parameter: The traffic demand is weighted precisely to offer few user or low traffic-demanding cells, at some point the chance of being illuminated. 
frame_dur=0.1;%[s] Each frame corresponds to 10 5G NR frames (15kHz SCS)=10*10ms=100ms.
frame=10; % Number of Illumination Slots. The total simulation time is fixed by the product between, frame*frame_dur: 100*0.1s=10s in this case. 
freq=20e9; % Central Operating Frequency: 20GHz -> Ka Band (DL)
h_sat=1100; % Starklink [550km] ; OneWeb [1100km] 
el_min=25; % Minimum Elevation Angle [º] 
Re=6371; % Earth Radius [km]
rings=10; % Number of cell rings to be evaluated within the FIXED SCENARIO
P_T=18; % Total RF power [W]
B_T=250/2; % Bandwidth per colour [MHz]
n_users=50; % Number of users within the satellite's Field of View (FoV).

% Trigonometric relationships for satellite's FoV determination:
betta_max=(asin((Re/(Re+h_sat))*cos(el_min*pi/180)))*180/pi;
gamma_max=90-betta_max-el_min;
h_prime_max=Re*(1-cos(gamma_max*pi/180));
A_footprint=2*pi*Re*h_prime_max;
D_footprint=((2*gamma_max)/360)*(2*pi*Re);

%FOM Parameter initialization:
RC=0; % Requested Capacity
SC=0;% Served Capaciaty (from requested)
UC=0; % Unserved Capacity
EC=0; % Extra Served Capacity
TTS=0; % Time To Serve

% Traffic Distribution - User generation:
traffic_model='random'; %'linear' %'hotspot'
[x_y,demand,type,g_rx,T_noise_rx]=Traffic_Distribution(traffic_model,n_users,D_footprint, freq);

% FoM Calculation:
cell_scenario_model='fixed'; %'variable'
[RC,SC,UC,EC,TTS]=BH_calculation(rings, P_T, B_T, n_users, cell_scenario_model, D_footprint, beams, colours, frame, frame_dur, TTL, freq, h_sat, el_min, x_y, demand, type, g_rx, T_noise_rx);

%Result Saving:
save(strcat('[',string(realization),']data_',string(beams),'beams_',string(colours),'colours',string(n_users),'users',string(P_T),'power'),"SC","UC","EC","RC","TTS")

% GIF Generation:
if display==1
    filename = convertCharsToStrings(strcat('[GIF]_',num2str(cell_scenario_model),'_cell_illumination_',num2str(rings),'rings.gif')); % Specify the output file name
    for i=2:frame+2
        figure (i)
           frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 2;
             imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.8);
        else
             imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.8);
        end
    end
    
    % Display Illumination Sequence:
    filename = convertCharsToStrings(strcat('[GIF]_MAP_',num2str(cell_scenario_model),'_cell_illumination_',num2str(rings),'rings.gif'));; % Specify the output file name
    for i=13:frame+10+2
        figure (i)
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 13;
             imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.8);
        else
             imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.8);
        end
    end
end

