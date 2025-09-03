function main
    % INPUTS
    
    % Paper vs simulation code input data: --------------------------------
    
    % paper has:
    % w_m = 0.7; n_m = 0.125; w_e = 3.5; n_e = 0.62
    % given inputs to Stewart et al. model are:
    % w_m = 3.5; n_m = 0.0273; w_e = 17.5; n_e = 0.0273
    % using the values from the model
    
    % max_growth is maximum growth rate for grass in a year
    % paper lists 0.125, while model takes 1.125
    % this is likely a difference in notation, since Creosote bush
    % similarly has 0.09 in paper and 1.09 in model input
    % however the model code appears to use 1.125*biomass as a limit for
    % propagules, so this one is used
    
    %  --------------------------------------------------------------------

    % take inputs from live script
    [Initial_Conditions, Field] = inputs;
    
    % constants for soil saturation approach
    water_saturation = 100;
    nitrogen_saturation = 50;
    
    % ---------- ! likewise for nitrogen
    d(1, 1:Initial_Conditions.T) = 1.5;
    
    % Grass properties ----------------------------------------------------
    
    % maximum biomass for grass in a cell 
    grass.b_max = 319;
    % add constants for plant behaviour
    grass.water_maintenance = 3.5;
    grass.nitrogen_maintenance = 0.0273;
    grass.water_efficiency = 17.5;
    grass.nitrogen_efficiency = 0.0273;
    % k is a total non-resource-related mortality rate, i.e. from diseases,
    % grazing, physical damage
    grass.k = 0.1;
    % f is failure rate for propagules to establish themselves
    grass.f = 0.05;
    % maximum growth per year
    grass.max_growth = 1.125;
    
    % ---------------------------------------------------------------------
    
    
    % DEFINE FUNCTIONS - temporarily put in the functions file
    
    
    % INITIALISE MODEL ----------------------------------------------------
    
    % setup rainfall, biomass record and deep layer resource record
    [Initial_Conditions, Field, rain] = initialise(Initial_Conditions, Field, grass);
    
    % MAIN ----------------------------------------------------------------
    biomass = Field.biomass_record(:, 1);
    deep_water = Field.deep_water_record(:, 1);
    deep_nitrogen = Field.deep_nitrogen_record(:, 1);
    for t = 1:Initial_Conditions.T
    
        % produce availability vectors for water and nitrogen
        water_availability = resource_availability(biomass, rain(t), Field.size, grass.b_max);
        nitrogen_availability = resource_availability(biomass, d(t), Field.size, grass.b_max);
    
        % generate matrix for local transport of resources and propagules
        sigma = generate_sig(biomass, Field.size, grass.b_max);
    
        % find final surface resource distributions using local transport
        % and availability
        surface_water = (rain(t) + sigma*water_availability);
        surface_nitrogen = d(t) + sigma*nitrogen_availability;
    
        % calculate WR and NR for use in biomass calculations
        water_res = (surface_water + deep_water - grass.water_maintenance.*biomass)./grass.water_efficiency;
        nitrogen_res = (surface_nitrogen + deep_nitrogen + grass.nitrogen_maintenance.*biomass)./grass.nitrogen_efficiency;
    
        % we construct a 3xN^2 matrix where each row is the resource 
        % remaining or maximum growth possible over the cells, by 
        % concatenating row vectors min returns the minimum for each column,
        % i.e. the 'I' value for each cell as a row vector, which we then 
        % change back to a column vector
        intermediate = min([water_res'; nitrogen_res'; grass.max_growth.*biomass'])';
    
        % any positive I values mean propagules are produced by a cell, 
        % otherwise we have zero biomass in a cell or insufficient resources 
        % and none are produced
        propagules = heaviside(intermediate).*intermediate;
    
        % calculate availability from this
        prop_availability = propagule_availability(biomass, propagules, Field.size, grass.b_max);
            
        % calculate insufficiency terms for the biomass updating equation
        % (this is not in the written form, but makes code more readable)
        water_insufficiency = (1-grass.k)*(grass.water_efficiency/grass.water_maintenance).*(1-heaviside(intermediate)).*heaviside(nitrogen_res-water_res).*intermediate;
        nitrogen_insufficiency = (1-grass.k)*(grass.nitrogen_efficiency/grass.nitrogen_maintenance).*(1-heaviside(intermediate)).*heaviside(water_res-nitrogen_res).*intermediate;
        
        % output some values for the central cell, used to investigate deep
        % soil accumulating too much
        % disp('At time t = ' + string(t+1) + ': B = ' + string(biomass(1250)) + '; W = ' + string(surface_water(1250)) + '; DW = ' + string(deep_water(1250)) + '; I = ' + string(intermediate(1250)))
    
        % update our soil resource stores 
        % if our cell is nearly empty, we remove resources to prevent
        % accumulation with:
        % (1-0.95.*heaviside(0.1*grass.b_max - biomass)).*
        deep_water = grass.water_efficiency.*(1-0.95.*heaviside(0.1*grass.b_max - biomass)).*(water_res - intermediate);
        deep_nitrogen = grass.nitrogen_efficiency.*(1-0.95.*heaviside(0.1*grass.b_max - biomass)).*(nitrogen_res - intermediate);
        
        % manually prevent too much resource from accumulating 
        % brute-force solution to resources accumulating too much in cells
        deep_water = water_saturation.*heaviside(deep_water - water_saturation) + heaviside(water_saturation - deep_water).*deep_water;
        deep_nitrogen = nitrogen_saturation.*heaviside(deep_nitrogen - nitrogen_saturation) + heaviside(nitrogen_saturation - deep_nitrogen).*deep_nitrogen;
    
        % update biomass
        biomass = (1-grass.k)*biomass + (1-grass.f)*propagules + (1-grass.f)*(sigma*prop_availability) + water_insufficiency + nitrogen_insufficiency;
    
        % if any biomass is above the maximum, reduce it to that
        % (note that floating point errors are protected against inside the
        % availability function, so we can set directly to maximum, rather 
        % than just below it)
        biomass = grass.b_max.*heaviside(biomass-grass.b_max) + heaviside(grass.b_max-biomass).*biomass;
    
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
    
    % OUTPUTS
    maps.grassmap = uint8([240 200 60; 230 200 60; 220 200 60; 210 200 60; 200 200 60; 180 200 60; 160 200 60; 140 200 60; 120 200 60]);
    maps.watermap = uint8([60 200 240; 50 180 220; 40 160 200; 30 140 180; 20 120 180; 10 100 160; 0 80 140]);
    maps.nitrogenmap = uint8([50 30 20; 80 40 18; 110 50 16; 140 60 14; 170 70 12; 200 80 10; 230 90 10]);
    
    % calculate an appropriate tick size for axes
    tick_size = (Field.size - mod(Field.size, 5))/5;
    if tick_size == 0
        tick_size = 1;
    end
    
    % find appropriate separation for ten timed outputs
    time_diff = (Initial_Conditions.T - mod(Initial_Conditions.T, 10))/10;
    if time_diff == 0
        time_diff = 1;
    end
    
    % PLOT GRAPHS ---------------------------------------------------------
    
    for t=1:time_diff:Initial_Conditions.T+1
        tick_graphs(Initial_Conditions, Field, maps, tick_size, t)
    end
    
    final_graphs(Initial_Conditions, Field, grass, rain, maps, tick_size, time_diff)
end