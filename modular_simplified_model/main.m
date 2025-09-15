function main

% ---------- NOTE: Paper vs simulation code input data ------------------ %

% paper has:
% w_m = 0.7; n_m = 0.125; w_e = 3.5; n_e = 0.62

% given inputs to Stewart et al. model are:
% w_m = 3.5; n_m = 0.0273; w_e = 17.5; n_e = 0.0273

% paper lists 0.125, while model takes 1.125 for max_growth
% this is likely a difference in notation, since Creosote bush similarly
% has 0.09 in paper and 1.09 in model input however the model code appears
% to use 1.125*biomass as a limit for propagules, so this one is used

% ---------- INITIALISE MODEL ------------------------------------------- %

% take inputs from live script
[Initial_Conditions, Field, Grass, Vectors, Graph_Options] = inputs;

% constants for soil saturation approach
water_saturation = 100;
nitrogen_saturation = 5;

% likewise for nitrogen
d(1, 1:Initial_Conditions.T) = 1.5;

% setup rainfall, biomass record and deep layer resource record
[Initial_Conditions, Field, Vectors, rain] = initialise(Initial_Conditions, Field, Grass, Vectors);

% ---------- MAIN ------------------------------------------------------- %

% set field records
biomass = Field.biomass_record(:, 1);
deep_water = Field.deep_water_record(:, 1);
deep_nitrogen = Field.deep_nitrogen_record(:, 1);

% create transport constants for each matrix
water_empty_transport_constant = create_transport_constants(Field.size, Vectors.Water.empty_transport);
water_grass_transport_constant = create_transport_constants(Field.size, Vectors.Water.grass_transport);
wind_empty_transport_constant = create_transport_constants(Field.size, Vectors.Wind.empty_transport);
wind_grass_transport_constant = create_transport_constants(Field.size, Vectors.Wind.grass_transport);

for t = 1:Initial_Conditions.T

    % produce availability vectors for water and nitrogen
    water_availability = resource_availability(biomass, rain(t), Field.size, Grass.b_max, Initial_Conditions.availability_periodicity);
    nitrogen_availability = resource_availability(biomass, d(t), Field.size, Grass.b_max, Initial_Conditions.availability_periodicity);

    % generate matrix for local transport of resources and propagules
    water_sigma = generate_sig(biomass, Grass.b_max, water_empty_transport_constant, water_grass_transport_constant);
    wind_sigma = generate_sig(biomass, Grass.b_max, wind_empty_transport_constant, wind_grass_transport_constant);

    % find final surface resource distributions using local transport
    % and availability
    surface_water = rain(t) + water_sigma*water_availability;
    surface_nitrogen = d(t) + water_sigma*nitrogen_availability*Vectors.Water.influence_index + wind_sigma*nitrogen_availability*Vectors.Wind.influence_index;

    % calculate WR and NR for use in biomass calculations
    water_res = (surface_water + deep_water - Grass.water_maintenance.*biomass)./Grass.water_efficiency;
    nitrogen_res = (surface_nitrogen + deep_nitrogen + Grass.nitrogen_maintenance.*biomass)./Grass.nitrogen_efficiency;

    % we construct a 3xN^2 matrix where each row is the resource 
    % remaining or maximum growth possible over the cells, by 
    % concatenating row vectors min returns the minimum for each column,
    % i.e. the 'I' value for each cell as a row vector, which we then 
    % change back to a column vector
    intermediate = min([water_res'; nitrogen_res'; Grass.max_growth.*biomass'])';

    % any positive I values mean propagules are produced by a cell, 
    % otherwise we have zero biomass in a cell or insufficient resources 
    % and none are produced
    propagules = heaviside(intermediate).*intermediate;

    % calculate availability from this
    prop_availability = propagule_availability(biomass, propagules, Field.size, Grass.b_max, Initial_Conditions.availability_periodicity);
        
    % calculate insufficiency terms for the biomass updating equation
    % (this is not in the written form, but makes code more readable)
    water_insufficiency = (1-Grass.k)*(Grass.water_efficiency/Grass.water_maintenance).*(1-heaviside(intermediate)).*heaviside(nitrogen_res-water_res).*intermediate;
    nitrogen_insufficiency = (1-Grass.k)*(Grass.nitrogen_efficiency/Grass.nitrogen_maintenance).*(1-heaviside(intermediate)).*heaviside(water_res-nitrogen_res).*intermediate;
    
    % output some values for the central cell, used to investigate deep
    % soil accumulating too much
    % disp('At time t = ' + string(t+1) + ': B = ' + string(biomass(1250)) + '; W = ' + string(surface_water(1250)) + '; DW = ' + string(deep_water(1250)) + '; I = ' + string(intermediate(1250)))

    % update our soil resource stores 
    % if our cell is nearly empty, we remove resources to prevent
    % accumulation with:
    % (1-0.95.*heaviside(0.1*Grass.b_max - biomass)).*
    deep_water = Grass.water_efficiency.*(1-0.95.*heaviside(0.1*Grass.b_max - biomass)).*(water_res - intermediate);
    deep_nitrogen = Grass.nitrogen_efficiency.*(1-0.95.*heaviside(0.1*Grass.b_max - biomass)).*(nitrogen_res - intermediate);
    
    % manually prevent too much resource from accumulating 
    % brute-force solution to resources accumulating too much in cells
    deep_water = water_saturation.*heaviside(deep_water - water_saturation) + heaviside(water_saturation - deep_water).*deep_water;
    deep_nitrogen = nitrogen_saturation.*heaviside(deep_nitrogen - nitrogen_saturation) + heaviside(nitrogen_saturation - deep_nitrogen).*deep_nitrogen;

    % update biomass
    biomass = (1-Grass.k)*biomass + (1-Grass.f)*propagules + (1-Grass.f)*(water_sigma*prop_availability*Vectors.Water.influence_index + wind_sigma*prop_availability*Vectors.Wind.influence_index) + water_insufficiency + nitrogen_insufficiency;

    % if any biomass is above the maximum, reduce it to that
    % (note that floating point errors are protected against inside the
    % availability function, so we can set directly to maximum, rather 
    % than just below it)
    biomass = Grass.b_max.*heaviside(biomass-Grass.b_max) + heaviside(Grass.b_max-biomass).*biomass;

    % ensure that biomass is non-negative
    % (seen a couple instances of a single cell spiralling to large 
    % negative values. may be caused by memory problem or a float error 
    % from biomass being near zero)
    biomass = heaviside(biomass).*biomass;

    % add values to the record array as a new column
    Field.biomass_record(:, t+1) = biomass;
    Field.deep_water_record(:, t+1) = deep_water;
    Field.deep_nitrogen_record(:, t+1) = deep_nitrogen;

end

% ---------- PLOT GRAPHS ------------------------------------------------ %

% colours for plots
maps.grassmap = uint8([240 200 60; 230 200 60; 220 200 60; 210 200 60; 200 200 60; 180 200 60; 160 200 60; 140 200 60; 120 200 60]);
maps.watermap = uint8([60 200 240; 50 180 220; 40 160 200; 30 140 180; 20 120 180; 10 100 160; 0 80 140]);
maps.nitrogenmap = uint8([50 30 20; 80 40 18; 110 50 16; 140 60 14; 170 70 12; 200 80 10; 230 90 10]);

% calculate an appropriate tick size for axes
tick_size = (Field.size - mod(Field.size, 5))/5;
if tick_size == 0; tick_size = 1; end

% find appropriate separation for ten timed outputs
time_diff = (Initial_Conditions.T - mod(Initial_Conditions.T, 10))/10;
if time_diff == 0; time_diff = 1; end

% tick graphs
if Graph_Options.plot_tick_graphs
    for t=1:time_diff:Initial_Conditions.T+1
        tick_graphs(Initial_Conditions, Field, maps, tick_size, t)
    end
end

% final graphs
final_graphs(Graph_Options, Initial_Conditions, Field, Grass, rain, maps, tick_size, time_diff)

end
