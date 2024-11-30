# Demand-based-Satellite-Beam-Hopping-BH-Algorithm

This code is © Samuel M. Zamacola, 2024, and it is made available under the GPL license enclosed with the software.
Over and above the legal restrictions imposed by this license, if you use this software for an academic publication then you are obliged to provide proper attribution to the paper that describes it:
+ S. M. Zamacola, R. M. Rodríguez-Osorio, and M. A. Salas-Natera, ‘Joint satellite platform and constellation sizing for instantaneous beam-hopping in 5G/6G Non-Terrestrial Networks’, Computer Networks, p. 110942, 2024. 
 
DEF.: "The demand-based algorithm is a deterministic algorithm where illumination scheduling (BH) and resource allocation is performed in response to users' traffic demand. It is dynamically adjusted according to the real-time needs of users at any given time. The overall flowchart is presented in Simulation_Flow.jpg. The diagram is divided in two parts: pre-computation modules to ease the assignment of resources and the calculation of the parameters of interest (PoI): UC/EC/TTS, and the BH algorithm to dynamically calculate illumination [Ill], power [P] and bandwidth [B] for each time slot. For more information, check the reference paper."
 
* INPUT: satellite altitude (h_sat), minimum elevation angle (el_min), number of rings within the satellite's FoV (rings), frequency (f), total ilumination slots (frame), slot duration (frame_dur), number of colours (colours), simultaneous beams (beams), ilumination weighting factor (TTL), number of users (n_users), distribution of users (traffic_model): random/linear/hotspot, cell scenario (cell_scenario_model): fixed/variable, total RF power (P_T), total bandwidth per colour (BW_T).

 
* OUTPUT (FoM): Unserved Capacity (UC), Extra Served Capacity (EC), Time To Serve (TTS).
 
Included Files:
 
+ User class (u.m): Each user is characterized by a given location, a type of station, a given traffic demand and
counts at UE level with some gain and noise characteristics. These are the static attributes of
the class, but there are other attributes that are going to be dynamically computed through
methods based on these prior attributes. Depending on the subsatellite point location and the
position of the satellite the following parameters are calculated: the slant range, the elevation
angle from the UE to the satellite, and the nadir angles, with respect to the satellite, and relative
to its cell center. Dynamic attributes, those that are going to be filled in each frame iteration and
therefore are time dependent, are used to store the intermediate and output variables of the
iteration. These are: assigned power and bandwidth resources for link budget computation, C/N
which is link budget’s output, and then served, pending and extra traffic counters.
 
+ Cell class (c.m): Cells are labelled by a given identifier or cell number, count with a given center and have a
radius that is going to be dependant of the total footprint area and the number of chosen rings.
Cell objects will include a position based computed attribute in which the user objects presented
above are stored. This is where the composition property comes into action, as cells count with
users located within their area, and the organisation of these is much simpler by incorporating
a list of user objects, if any, in each of the cells. Additionally, the nadir angle from the cell centre
and the list of interfering cells is computed based on the initial scheme. During each frame
iteration the aggregated traffic requirement is calculated based on the user list within the given
cell. Based on the cells that are selected for illumination, these are going to be classified in each
time instant as active or inactive. If that is the case, the aggregated power and bandwidth
resources assigned at cell level are computed. In case the cell is active, the assigned colour is
also stored as a dynamic attribute.
 
+ main.m: the executable file where BH algorithm is executed. Input parameters are defined in the file and both pre-computation (traffic_model->Traffic_Distribution, cell_scenario_model->Cell_Scenario), and BH computation (BH_Calculation) funtions are called through it.
 
+ Traffic_Distribution.m: based on the selected user distribution type, the generation of users is performed in the file, by defining UE related specifications.
 
+ Cell_Scenario.m: based on the number of rings that are intended to be allocated within the satellite's FoV, cells are generated in the file.

+ Non_Fixed_Beam_Aggrupation.m: used for the generation of cells in a variable manner (cell_scenario_model) following the Bron-Kerbosch maximal independent set and maximal clique algorithm. IMPORTANT! In addition to the fucntion, for variable cell generation, 'BK_MaxClique.m' and 'BK_MaxIS.m' auxiliary functions must be added into the execution folder. Find them in: https://es.mathworks.com/matlabcentral/fileexchange/24591-bron-kerbosch-maximal-independent-set-and-maximal-clique-algorithms

+ BH_Calculation.m: execution of the BH algorithm. Once these pre-computations are preformed (users and cells), for each time slot, cell illumination [Ill], resource allocation ([B] and [P]) and link budget calculations are performed, determining the FoM: EC, UC, TTS. To consider EXTRA losses rather than FSL ones, within the link budget equation, uncomment lines 217-252. The following files are required to be downloaded and added into the execution folder: 'maps.mat','p836.mat','p837.mat' and 'p840.mat'. Find them in: https://es.mathworks.com/help/satcom/ref/p618propagationlosses.html.
 
NOTE: For a fast execution of the algorithm set display=0 in main.m. No figures/GIFs will be displayed. Fixed and Variablexample GIFs: [GIF]_Illumination_xx_rings.gif and [GIF]_Illumination_MAP_xx_rings.gif 

## Contact
For questions, issues, or contributions, please contact Samuel M. Zamacola at samuel.martinez@upm.es.
