function [RC,SC,UC,EC,TTS]=BH_calculation(rings, P_T, B_T, n_users, cell_scenario_model, D_footprint, beams, colours, frame, frame_dur,  TTL, freq, h_sat, el_min, x_y, demand, type, g_rx, T_noise_rx)

global display;

RC=0;
SC=0;
UC=0;
EC=0;
TTS=0;


% CQI Table 1: ACM
SE_t1=[0.1523, 0.2344, 0.3770, 0.6016, 0.8770, 1.1758, 1.4766, 1.9141, 2.4063, 2.7305, 3.3223, 3.9023, 4.5234, 5.1152, 5.5574];
% CQI Table 2: ACM
SE_t2=[0.1523, 0.3770, 0.8770, 1.4766, 1.9141, 2.4063, 2.7305, 3.3223, 3.9023, 4.5234, 5.1152, 5.5547, 6.2266, 6.9141, 7.4063];

if strcmp(cell_scenario_model,'fixed')
    %% FIXED Cell Scenario
    if display==1
        % PAINT
        % Subsatellite Point Selection:
        figure (1)
        axesm robinson
        gridm on
        geoshow('landareas.shp')
        title('Select Subsatellite Point:')
        % [lat_sp, lon_sp] = inputm(1);
        lat_sp=50; lon_sp=-30;
        geoshow(lat_sp, lon_sp, 'DisplayType', 'Point', 'Marker', 'o', 'MarkerEdgeColor', 'black','MarkerFaceColor','black','MarkerSize',5)
        pause (1)
    end
    
    % Single patch antenna radiation pattern cosq(theta) approximation:
    % syms q
    % q=vpasolve(cosd(theta/2)^q==1/sqrt(2));
    % q=double(q);
    q=1.3;
    
    np=1001;
    theta_scanning=linspace(-pi/2,pi/2,np);
    % Eo=1;
    % Ep=Eo*(cos(theta_scanning)).^q;
    % figure 
    % plot(rad2deg(theta_scanning),20*log10(Ep)); 
    % xlim([-90 90]);
    % ylim([-10 0]);
    % title('Element radiation pattern')
    % xlabel('theta (°)')
    % ylabel('dB')
    
    
    % Cell Definition:
    [l,h,theta,D_footprint,number_cells,centers, interfering]=Cell_Scenario(h_sat,el_min,rings);
    
    if display==1
        % PAINT
        % Satellite Footprint:
        figure (2)
        pgon = nsidedpoly(1000, 'Center',  [0,0], 'Radius', D_footprint/2);
        plot(pgon, 'FaceColor', 'y');
        axis equal
        hold on
    end
    
    % Cell Footprint :
    c_scenario=[];
    for i=1:number_cells
        c_scenario=[c_scenario c(i,centers(i,:),l,nonzeros(interfering(i,:)))]; %constructor
        c_scenario(i).compute_betta_to_sat(h_sat); % Compute betta to the cell center so that the gain loss due to beam scanning (moving ftom boresight) can be then accounted: non-ideal isotropic behavior of the embedded element gain.
       if display==1
           % PAINT
           c_scenario(i).draw(1);
           hold on
       end
    end
    if display==1
        % PAINT
        title(strcat('Cell Scenario:  ',num2str(rings),' rings'))
    end
    
    
    % Cell Assignment:
    for i=1:n_users
        % Euclidean distance to cell centers to determine the cell:
        [dist_min,cell_num]=min(sqrt((x_y(i,1)-centers(:,1)).^2+(x_y(i,2)-centers(:,2)).^2)); 
        % Assign user to corresponding cell:
        c_scenario(cell_num).adduser(u([x_y(i,1),x_y(i,2)],type,demand(i),g_rx,T_noise_rx));  %constructor
        c_scenario(cell_num).users(length(c_scenario(cell_num).users)).compute_distance_elevation_betta_to_sat(h_sat); % Compute distance, elevation and betta angle to the satellite of the the last added user
        c_scenario(cell_num).users(length(c_scenario(cell_num).users)).compute_betta_to_cell_center(h_sat,c_scenario(cell_num)); % Compute the betta angle with respect to the cell center so that the gain loss can be accounted by the fact of not being at the cell center
    end

elseif strcmp(cell_scenario_model,'variable')
    %% VARIABLE Cell Scenario
    if display==1
    % PAINT
    % Subsatellite Point Selection:
    figure (1)
    axesm robinson
    gridm on
    geoshow('landareas.shp')
    title('Select Subsatellite Point:')
    % [lat_sp, lon_sp] = inputm(1);
    lat_sp=50; lon_sp=-30;
    geoshow(lat_sp, lon_sp, 'DisplayType', 'Point', 'Marker', 'o', 'MarkerEdgeColor', 'black','MarkerFaceColor','black','MarkerSize',5)
    pause (1)
    end
    
    % Cosq(theta) approximation:
    % syms q
    % q=vpasolve(cosd(theta/2)^q==1/sqrt(2));
    % q=double(q);
    q=1.3;
    
    np=1001;
    theta_scanning=linspace(-pi/2,pi/2,np);
    % Eo=1;
    % Ep=Eo*(cos(theta_scanning)).^q;
    % figure 
    % plot(rad2deg(theta_scanning),20*log10(Ep)); 
    % xlim([-90 90]);
    % ylim([-10 0]);
    % title('Element radiation pattern')
    % xlabel('theta (°)')
    % ylabel('dB')
    
    users=[];
    for i=1:n_users
    users=[users u([x_y(i,1),x_y(i,2)],type,demand(i),g_rx,T_noise_rx)];
    users(i).compute_distance_elevation_betta_to_sat(h_sat);
    end
    
    % Calculate FoV and l values based on rings:
    [l,h,theta]=Cell_Scenario(h_sat,el_min,rings);
    FoV=theta; 
    
    if display==1
        % PAINT
        % Satellite Footprint:
        figure (2)
        pgon = nsidedpoly(1000, 'Center',  [0,0], 'Radius', D_footprint/2);
        plot(pgon, 'FaceColor', 'y');
        axis equal
        hold on
    end
    
    adj_matrix=zeros(n_users,n_users); % Adjacent matrix: 1 if users can be allocated into the same beam (alpha(i,j)<FoV), if not 0
    
    % Distance between users:
    dist=zeros(n_users,n_users);
    
    for i=1:n_users
        for j=1:n_users
            dist=sqrt((users(i).location(1)-users(j).location(1))^2+(users(i).location(2)-users(j).location(2))^2);
            if dist<=2*l
                adj_matrix(i,j)=1;
            end
        end
    end 
    
    % Non Fixed Beam - User Aggrupation Algorithm:
    algorithm_repetition=1; % Increase the value to avoid local optima!
    S=Non_Fixed_Beam_Aggrupation(n_users, adj_matrix, algorithm_repetition);
    
    % Cell definition by calculating the centroid of each user aggrupation:
    number_cells=length(S);
    % Cell Footprint:
    c_scenario=[];
    centers=zeros(number_cells,2);
    for i=1:number_cells
        % Centroid calculation;
        for user_idx=1:length(S{i})
            centers(i,1)=centers(i,1)+users(S{i}(user_idx)).location(1);
            centers(i,2)=centers(i,2)+users(S{i}(user_idx)).location(2);
        end
        centers(i,1)=centers(i,1)/length(S{i});
        centers(i,2)=centers(i,2)/length(S{i});
       
        c_scenario=[c_scenario c(i,centers(i,:),l,[])]; %constructor: for this case assume interfering []
        c_scenario(i).compute_betta_to_sat(h_sat); % Compute betta to the cell center so that the gain loss due to beam scanning (moving ftom boresight) can be then accounted: non-ideal isotropic behavior of the embedded element gain.
        
        % Assign user to corresponding cell:
        for user_idx=1:length(S{i})
            c_scenario(i).adduser(users(S{i}(user_idx)));  %constructor
            c_scenario(i).users(user_idx).compute_distance_elevation_betta_to_sat(h_sat); % Compute distance, elevation and betta angle to the satellite of the the last added user
            c_scenario(i).users(user_idx).compute_betta_to_cell_center(h_sat,c_scenario(i)); % Compute the betta angle with respect to the cell center so that the gain loss can be accounted by the fact of not being at the cell center
        end
        if display==1
           % PAINT
           hold on
           c_scenario(i).draw(1);
        end
    end
   
        % Interfering Cells Computation Variable Case (do not consider the output from Cell_Scenario):
        d_interfering=4*h;
        interfering=zeros(number_cells,1);
        for i=1:number_cells
            cont=0;
            for j=1:number_cells
                if i~=j
                    if sqrt((c_scenario(i).center(1)-c_scenario(j).center(1))^2+(c_scenario(i).center(2)-c_scenario(j).center(2))^2)<d_interfering
                        cont=cont+1;
                        interfering(i,cont)=j;
                    end
                end
            end
        end
    
    
        for i=1:number_cells
            c_scenario(i).interfering=nonzeros(interfering(i,:));
        end
else 
    fprintf('No correct cell scenario option selected!\n')

end
%% Traffic Flow:
cell_illumination_percentage=zeros(1,length(c_scenario));
P_t=10*log10(P_T);
G_t=10*log10(0.65*48360/theta^2);
%P_t=Power_Reduction-G_t-10*log10(beams); %dBW %10*log10(80); 
k=-228.6012; % Boltzmann Constant [dBW/(HzK^-1)]
T_ant=30; % Antenna Temperature [K]
alpha=0.1; % Rolloffactor


for t=1:frame  % Illumination switching performed each frame: as frame_dur=0.1, for 10 5G frames constant illumination.

    % Traffic Pending Adjustment based on frame length: from second to second add the total demand/s to the pending:
    if mod(t-1,1/frame_dur)==0
        for i=1:length(c_scenario)
            for j=1:length(c_scenario(i).users)
                c_scenario(i).users(j).traffic_pending=c_scenario(i).users(j).traffic_pending+c_scenario(i).users(j).traffic_demand;
            end
        end
    end

    % Agregated_Traffic:
    cell_traffic=zeros(1,length(c_scenario));
    for c_idx=1:length(c_scenario)
        c_scenario(c_idx).aggregatetraffic(); % Compute the aggregated cell traffic through the method.
        cell_traffic(c_idx)=c_scenario(c_idx).aggregated_traffic(t); % Take the property that has been generated by the previous method.
        cell_illumination_percentage(c_idx)=sum(c_scenario(c_idx).active)/length(c_scenario(c_idx).active); % Compute the percentage of previous tramas in which the respective (c) cell has been illuminated.
        c_scenario(c_idx).active(t)=0;
        c_scenario(c_idx).colour(t)=0;
        c_scenario(c_idx).P(t)=0;
        c_scenario(c_idx).BW(t)=0;
        for u_idx=1:length(c_scenario(c_idx).users)
            c_scenario(c_idx).users(u_idx).P(t)=0;
            c_scenario(c_idx).users(u_idx).BW(t)=0;
            c_scenario(c_idx).users(u_idx).C_N(t)=0;
            c_scenario(c_idx).users(u_idx).traffic_served(t)=0;
        end
        if (t>=TTL && sum(c_scenario(c_idx).active(t-TTL+1:t))==0)
            cell_traffic(c_idx)=cell_traffic(c_idx)*(1+(1-cell_illumination_percentage(c_idx))); % Weighted traffic aggregation: if it has not been previously selected increase its weight.
        end     
    end

    if display==1
        % PAINT
        % Map plot per slot:
        figure (t+2+10) % MAP FIGURE
        geoshow('landareas.shp')
        % Satellite Footprint:
        hold on
        pgon = nsidedpoly(1000, 'Center',  [0,0], 'Radius', D_footprint/2);
        Footprint_Latitude=lat_sp+(pgon.Vertices(:,2))/110.574; % approximated conversion from km diff. in x and y to degrees in lat. and long.
        Footprint_Longitude=lon_sp+(pgon.Vertices(:,1))./(111.32.*cosd(Footprint_Latitude));
        plot(polyshape(Footprint_Longitude',Footprint_Latitude'),'FaceColor','yellow','EdgeColor', 'black')
    end
        
    % Illumination Sequence (MAXIMUM ORDERING + WEIGHTED SELECTION):
    [max2min_ni,idx_max2min_ni]=sort(cell_traffic,'descend'); % It could happen that the values ordered in max2min do not match with the aggregated pending traffic, as it can be ponderated/weighted depending on the TTL.

    % Co-Channel Interference (colour) restriction:
    max2min=zeros(1,length(max2min_ni));
    idx_max2min=zeros(1,length(idx_max2min_ni));
    colour=[];
    beam_idx=0;
    prohibited_cell=0;
    updated_beams=beams;
    while (beam_idx<updated_beams) && (beam_idx<number_cells)
        beam_idx=beam_idx+1;
            if beam_idx<=colours % CO-CHANNEL INTERFERENCE! No colour re-use:
                max2min(beam_idx)=max2min_ni(beam_idx);
                idx_max2min(beam_idx)=idx_max2min_ni(beam_idx);
                colour(beam_idx)=beam_idx;
            else % CO-CHANNEL INTERFERENCE! We need to re-use colour:
                free_colours=[];
                prohibited_colours=[];
                for colour_idx=1:length(colour) % Check past ones.
                    if ~(any(c_scenario(idx_max2min_ni(beam_idx)).interfering == idx_max2min(colour_idx))) && ~any(prohibited_colours == colour(colour_idx)) % If the next cell that we are willing to illuminate is not limited in the colour that can be chosen, or in other words, if the cell does not contain in his interfering list the cell number whose colour has already been assigned; AND, if the colour is not in the prohibited colour list by previous assigments:
                        free_colours=[free_colours colour(colour_idx)];
                    else  % Comparisons are done 1 to 1 so if a cell has been considered interfering one, its colour should be marked as prohibited so that it can not be chosen: 
                        prohibited_colours=[prohibited_colours colour(colour_idx)]; %  It should be removed from the free colour list. 
                    end
                end
                for idx_pc=1:length(prohibited_colours) % Take out from free_colours any prohibited_colours that could exist.
                    free_colours=free_colours(free_colours~=prohibited_colours(idx_pc));
                end
                if ~isempty(free_colours) % If there is a free colour:
                   % The assigment between the free colours should be to select the one that is less frequently chossen, the minimum repeated one:
                   unique_colours=unique(free_colours); % unique colours present in the free_colours vector.
                   repetitiveness_colours= histc(free_colours,unique_colours); % repetitiveness in the vector of the unique values of that same vector.
                   [rc_min,rc_idx] = min(repetitiveness_colours); % Minimum ordenation.
                   colour(beam_idx-prohibited_cell) = unique_colours(rc_idx); % Assign the minimum one.
                   max2min(beam_idx-prohibited_cell)=max2min_ni(beam_idx);
                   idx_max2min(beam_idx-prohibited_cell)=idx_max2min_ni(beam_idx);
                elseif isempty(free_colours) && beam_idx<number_cells % It is the case of a prohibited cell! There is no extra colour available, let's give the chance to the next maximum aggregated cell. Searching can continue as there are more cells to illuminate (if the condition is satisfied).
                    updated_beams=updated_beams+1; % Extra beam possible as the previous cell has not been illuminated.
                    prohibited_cell=prohibited_cell+1;
                end
            end
    end
    % Restriction if the number of beams is greater that the cells that in each ring scenario exist:
    if number_cells<beams
        beams=number_cells;
    end

    % Resource Allocation per Cell:
    for beam_idx=1:beams % Select the already ordered and interference limitted beam values.
        if max2min(beam_idx)==0 % If the aggregated traffic if 0 then:
            break % No need to keep illuminating!
        else % Perform the calculations:
            c_scenario(idx_max2min(beam_idx)).active(t)=1;
            c_scenario(idx_max2min(beam_idx)).colour(t)=colour(beam_idx);
            % Power and BW distribution at cell level:
            c_scenario(idx_max2min(beam_idx)).P(t)=P_t+10*log10(max2min(beam_idx)/sum(max2min(1:beams))); % [dBW] not uniform distribution, total power is distributed in each beam for the whole trama based on aggregated traffic demand <- INCORPORAR CRITERIO G/T !!!
            c_scenario(idx_max2min(beam_idx)).BW(t)=10*log10(B_T*10^6); %[dBHz]
            % For each user in the selected colored cell and time slot, compute:
            for u_idx=1:length(c_scenario(idx_max2min(beam_idx)).users)
                if (c_scenario(idx_max2min(beam_idx)).users(u_idx).traffic_pending~=0) % If ==0: no need to serve this user of the cell, the illumination is to serve other demanding users of the cell.
                    % Link Budget:
                    % C/N=Prx-Pnoise=(EIRPsat-Lfsl-Lextra+Grx)-(k+Teq+B)  |dB FOR EACH FRAME: (10SLOTS)
                    % EIRP_tx=P_tx+G_tx:
                    P_tx=c_scenario(idx_max2min(beam_idx)).P(t)+10*log10(c_scenario(idx_max2min(beam_idx)).users(u_idx).traffic_pending/c_scenario(idx_max2min(beam_idx)).aggregated_traffic(t)); % The P of the cell, which has already been distributed from the total one (P_t) based on the aggregated traffic, is now divided based on the traffic required (pending) by the users in the cell to whom resources will be allocated.
                    G_tx=G_t-12*(c_scenario(idx_max2min(beam_idx)).users(u_idx).betta_to_center/theta)^2; % Gain loss approximation when moving away from the cell center.
                    [minScanning,closestScanningIndex]=min(abs(theta_scanning(1:((np-1)/2+1))+deg2rad(c_scenario(idx_max2min(beam_idx)).betta_to_sat))); % Gain loss due to scanning
                    G_tx=G_tx+20*log10(cos(abs(theta_scanning(closestScanningIndex)))^q);
                    % FSL:
                    Lfsl=20*log10(freq)+20*log10(c_scenario(idx_max2min(beam_idx)).users(u_idx).distace_to_sat*10^3)-147.55;

%                   EXTRA LOSSES: https://es.mathworks.com/help/satcom/ref/p618propagationlosses.html
%                   DATA:
%                            % "maps.m": lat/long grid map:
%                                 % R001: rainfall rate exceeded for 0.01% of average year. [mm/h]
%                                 % p839 -> h0: mean anual 0ºC isotherm height AMSL [km]
%                                 % p453->NWET_Annual_50: wet term surface refractivity exceeded for 50% of the year
%                                 % p1511->alt_file: topographic height AMSL
%                                 % P1510->temp_file: mean surface temperature at 2m above the surface of the Earth
%           
%                             % “p836.m”: Rho: surface water vapour density [g/m3], Vsch: water vapour scale height [km], Vt: total water vapour content [kg/m2]
%                             % “p837.m”: Rp: rainfall rate exceeded for the desired anual probability p of excedance 
%                             % "p840.mf: Lred: columnar content of liquid water reducef to a temperature of 273.15K [kg/m2 or mm] at a probability p (clouds and fog)
% 
%                     cfg=p618Config;
%                     cfg.Frequency=freq;
%                     cfg.ElevationAngle=c_scenario(idx_max2min(beam_idx)).users(u).elevation_to_sat
%
%                     x
%                     cfg.Latitude=lat_sp+(c_scenario(idx_max2min(beam_idx)).users(u).location(2))/110.574; %approximated conversion from km diff. in x and y to degrees in lat. and long.
%                     cfg.Longitude=lon_sp+(c_scenario(idx_max2min(beam_idx)).users(u).location(1))/(111.32*cosd(cfg.Latitude));
%                     if display==1
%                       figure (t+2+10) %MAP FIGURE
%                       geoshow(cfg.Latitude, cfg.Longitude, 'DisplayType', 'Point', 'Marker', 'o','MarkerEdgeColor',char(c_scenario(idx_max2min(beam_idx)).colour_colours(c_scenario(idx_max2min(beam_idx)).colour(t))) ,'MarkerFaceColor',char(c_scenario(idx_max2min(beam_idx)).colour_colours(c_scenario(idx_max2min(beam_idx)).colour(t))),'MarkerSize',4)          
%                     end
%                     % Average annual time percentage of excess for...  
%                     % For an availability of 99.9%, the excess is of 0.1. 
%                     cfg.GasAnnualExceedance = 0.1; % [0.1, 99] Ag!
%                     cfg.CloudAnnualExceedance = 0.1; % [0.1, 99] Ac!
%                     cfg.RainAnnualExceedance = 0.1; % [0.001, 5] Ar!
%                     cfg.ScintillationAnnualExceedance=0.1; % [0.01, 50] As!
%                     cfg.TotalAnnualExceedance = 0.1; % [0.001, 50] At! -> For systems operating at frequencies above about 18 GHz, and especially those operating with low elevation angles and/or margins, the effect of multiple sources of simultaneously occurring atmospheric attenuation must be considered: AT(p)=AG(p)+sqrt((AR(p)+AC(p))^2+AS(p)^2)
%                     cfg.AntennaDiameter=(70*3*10^8)/(theta*cfg.Frequency);
%                     cfg.AntennaEfficiency=0.65;
%         
%                     Extra_losses=p618PropagationLosses(cfg); % See LINK to include LOCAL DATA arguments. If the local data is not available as an input, the function uses the digital maps provided in ITU-R to estimate each of the parameters: StationHeight, Temperature, Pressure, WaterVaporDensity, IntegratedWaterVaporContent, TotalColumnarContent, RainRate, WetSurfaceRefractivity, MedianRadiatingTemperature.
%                     % Extra_losses.Ag; Extra_losses.Ac; Extra_losses.Ar; Extra_losses.As; Extra_losses.At;
                    
                    % IF NO EXTRA LOSSES:
                    Extra_losses.At=0;
        
                    % Gain RX:
                    G_rx=c_scenario(idx_max2min(beam_idx)).users(u_idx).gain; %dB
                    % T_eq RX:
                    T_eq=10*log10(T_ant+c_scenario(idx_max2min(beam_idx)).users(u_idx).T_noise+280*(1-10^(-Extra_losses.At/10))); %[dBK] -> Contributors: captured by antenna + rx equivalent + attenuation due to extra losses
                    % BW RX:
                    BW_rx=c_scenario(idx_max2min(beam_idx)).BW(t)+10*log10(c_scenario(idx_max2min(beam_idx)).users(u_idx).traffic_pending/c_scenario(idx_max2min(beam_idx)).aggregated_traffic(t)); % As with the cell power, the cell BW divided based on the traffic demand of the users in the cell to whom resources will be allocated.
        
                    C_N=(G_tx+P_tx-Lfsl-Extra_losses.At+G_rx)-(k+T_eq+BW_rx);
                    
        
                    % Throughput Calculation (C/N -> bits/s): Adaptative ModulatioN and Coding Technology in 5G.
        
        %             % Table 1: For Block Error Rate of 10^-1 and AWWGN channel:
        %             if C_N<0
        %                 CQI=0.5633*C_N+4.9040;
        %             elseif C_N>=0 && C_N<7
        %                 CQI=0.5639*C_N+4.7443;
        %             elseif C_N>=7 && C_N<14
        %                 CQI=0.5493*C_N+4.5880;
        %             else
        %                 CQI=0.5319*C_N+4.5624;
        %                 if CQI>15 %not to exceed the dimensions of the CQI array, 15 is the maximum.
        %                     CQI=15;
        %                 end
        %             end
        
                    % Table 2:
                    if C_N<4
                        CQI=0.2789*C_N+2.9142;
                    elseif C_N>=4 && C_N<14
                        CQI=0.5246*C_N+1.8749;
                    else
                        CQI=0.5369*C_N+1.4946;
                        if CQI>15 %not to exceed the dimensions of the CQI array, 15 is the maximum.
                            CQI=15;
                        end
                    end
                    
                    if CQI>1
                        SE=SE_t2(floor(CQI));
                    else %  For really low values of C_N the link can not be ensured, as the required quality (CQI) for basic MODCODs can no be achieved. Therefore, SE=0.
                        SE=0;
                    end
                    Rs=(10^(BW_rx/10)/(10^6))/(1+alpha); % Symbol Rate [Msimbol/s] 
                    Throughput=Rs*SE; % Throughput [Msimbol/s*bit/simbol]=[Mbps]   
       
                    % Value Storing:
                    c_scenario(idx_max2min(beam_idx)).users(u_idx).P(t)=P_tx;
                    c_scenario(idx_max2min(beam_idx)).users(u_idx).BW(t)=BW_rx;
                    c_scenario(idx_max2min(beam_idx)).users(u_idx).C_N(t)=C_N;
                    c_scenario(idx_max2min(beam_idx)).users(u_idx).traffic_served(t)=Throughput;
                    c_scenario(idx_max2min(beam_idx)).users(u_idx).traffic_pending=c_scenario(idx_max2min(beam_idx)).users(u_idx).traffic_pending-Throughput/(10);
                    if c_scenario(idx_max2min(beam_idx)).users(u_idx).traffic_pending<0
                        EC=EC+(-c_scenario(idx_max2min(beam_idx)).users(u_idx).traffic_pending);
                        c_scenario(idx_max2min(beam_idx)).users(u_idx).traffic_pending=0;
                    end
                end
            end


        end 
    end
    
    if display==1
        % PAINT
        % Plot the cells (colours) that are going to be illuminated in this FRAME:
        figure (t+2+10) %MAP FIGURE
        title(append('Cell Scenario (t=',num2str(t),')'))
        figure (t+2) %SCHEMATIC FIGURE
        title(append('Cell Scenario (t=',num2str(t),')'))
        for i=1:number_cells
            figure (t+2) %SCHEMATIC FIGURE
            hold on
            c_scenario(i).draw(t);
            figure (t+2+10) %MAP FIGURE
            hold on
            c_scenario(i).draw_latlong(lat_sp,lon_sp,t);
            xlim([-80,20])
            ylim([0,90])
        end
        hold off
    end


end



%% FOM indicators: RC, SC, UC, EC

% Capacity:
for c_idx=1:length(c_scenario)
    for u_idx=1:length(c_scenario(c_idx).users)
        RC=RC+((frame/10)*c_scenario(c_idx).users(u_idx).traffic_demand);
        SC=SC+((frame/10)*c_scenario(c_idx).users(u_idx).traffic_demand-c_scenario(c_idx).users(u_idx).traffic_pending);
        UC=UC+c_scenario(c_idx).users(u_idx).traffic_pending;
    end
end

% Time:
TTS_cell=zeros(1,length(c_scenario));
for c_idx=1:length(c_scenario)
    act_cont=0;
    for idx=1:length(c_scenario(c_idx).active)
        if c_scenario(c_idx).active(idx)==1
            TTS_cell(c_idx)=TTS_cell(c_idx)+act_cont;
            act_cont=0;
        elseif c_scenario(c_idx).aggregated_traffic(idx)>0 % A cell will be waiting for illumination only if there is traffic to be served, if not, not fair to account the non serving, no users needing it!
            act_cont=act_cont+1;
 
        end
    end
    if c_scenario(c_idx).active(idx)==0 && c_scenario(c_idx).aggregated_traffic(idx)>0 % Consider last 0s as still more traffic to serve. It's best case indeed as it is not assured that in the next frame: frame_total+1 it is going to be served
            TTS_cell(c_idx)=(TTS_cell(c_idx)+act_cont)/(sum(c_scenario(c_idx).active)+1);
    elseif sum(c_scenario(c_idx).aggregated_traffic)>0 % Avoid NaN: the case where it's not last index being 1 but the one with no aggragated traffic ever.
        TTS_cell(c_idx)=(TTS_cell(c_idx))/(sum(c_scenario(c_idx).active));
    end
    % User ponderation:
    TTS_cell(c_idx)=TTS_cell(c_idx)*length(c_scenario(c_idx).users);
end
 
TTS=sum(TTS_cell)/n_users;
TTS=TTS*frame_dur; %TTS in absolute units: [s]

end