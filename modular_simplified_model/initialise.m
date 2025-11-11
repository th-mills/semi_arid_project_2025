function [Initial_Conditions, Field, Vectors, rain] = initialise(Initial_Conditions, Field, Grass, Vectors)

    % generate initial biomass distribution
    switch Initial_Conditions.initial_biomass_type

        case 0 % uniform random
            biomass = Grass.b_max*0.3*rand(Field.size^2,1);
    
        case 1 % stripes 
            biomass = 0.35*Grass.b_max*rand(Field.size,Field.size);
            for i = 1:Field.size
                biomass(i,:) = 0.5*(1+sin(16*pi*i/Field.size))*biomass(i,:);
            end
            biomass = reshape(biomass,Field.size^2,1);
    
        case 2 % spots
            biomass = 0.6*Grass.b_max*rand(Field.size,Field.size);
            for i = 1:Field.size
                biomass(i,:) = 0.5*(1+sin(10*pi*i/Field.size))*biomass(i,:);
            end
            for j = 1:Field.size
                biomass(:,j) = 0.5*(1+sin(10*pi*j/Field.size))*biomass(:,j);
            end
            biomass = reshape(biomass,Field.size^2, 1);

        case 3
            % central spot
            biomass = zeros(Field.size);
            biomass(0.4*Field.size:0.6*Field.size, 0.4*Field.size:0.6*Field.size) = 0.5*Grass.b_max*rand(0.2*Field.size+1);
            biomass = reshape(biomass, Field.size^2, 1);

        case 4
            % biomass around sides
            biomass = 0.3*Grass.b_max*rand(Field.size);
            biomass(0.2*Field.size:0.8*Field.size, 0.2*Field.size:0.8*Field.size) = 0;
            biomass = reshape(biomass, Field.size^2, 1);

        case 5
            % uniform random using random100.dat
            noise = readmatrix('random100.dat');
            noise = reshape(noise', [], 1);
            noise = noise(1:Field.size^2);
            biomass = 60*ones(Field.size^2,1).*noise();

        % case 3 % small patch in centre
        %     biomass = zeros(Field.size,Field.size);
        %     for i = Field.size/2-5:Field.size/2+5
        %         disp(i)
        %         for j = Field.size/2-5:Field.size/2+5
        %             biomass(i,j) = Grass.b_max*0.3;%*rand();
        %         end
        %     end
        %     biomass = reshape(biomass,Field.size^2, 1);
        % 
        % case 4 % small hole in centre
        %     biomass = 0.3*rand(Field.size,Field.size);
        %     for i = Field.size/2-5:Field.size/2+5
        %         disp(i)
        %         for j = Field.size/2-5:Field.size/2+5
        %             biomass(i,j) = 0;
        %         end
        %     end
        %     biomass = reshape(biomass,Field.size^2, 1);
    
    end
    
    % generate initial soil resource distribution
    % used 25 and 0.5 from initial values given in Stewart et al. 
    deep_water = zeros(Field.size^2, 1);
    deep_water(1:Field.size^2, 1) = 25;
    deep_nitrogen = zeros(Field.size^2, 1);
    deep_nitrogen(1:Field.size^2, 1) = 0.5;
    
    % rainfall uses an average of 243 mm/year by default
    rain = zeros(Initial_Conditions.T, 1);
    rain(1:Initial_Conditions.T, 1) = 243;
    Initial_Conditions.start_year = 0;

    switch Initial_Conditions.rain_type

        case 1 % historical rainfall
            % get the rain data as a matrix, with years in first column and
            % rainfall in second
            raindata = readmatrix("raindat.dat");
            % change our start year to that of the historical data
            Initial_Conditions.start_year = raindata(1, 1);
            if Initial_Conditions.T > size(raindata, 1)
                % if our simulation goes beyond the length of the record, change
                % only those up to the time-length of the record
                rain(1:size(raindata, 1)) = raindata(:, 2);
            else
                rain(1:Initial_Conditions.T) = raindata(1:Initial_Conditions.T, 2);
            end
        
        case 2 % repeated historical rainfall
            raindata = readmatrix("raindat.dat");
            % change our start year to that of the historical data
            Initial_Conditions.start_year = raindata(1, 1);
            len = size(raindata, 1);
            % get quotient and remainder of Initial_Conditions.T / raindata length
            remainder = mod(Initial_Conditions.T, len);
            quotient = (Initial_Conditions.T-remainder)/len;
            % add repeated historical record
            for n=0:(quotient-1)
                rain(n*len+1:(n+1)*len) = raindata(:, 2);
            end
            % add final section which is less than the full raindata length
            rain(quotient*len+1:quotient*len+remainder) = raindata(1:remainder, 2); 
    end
    
    % manually add droughts (fun to play with!)
    
    % 0.1 gives sparse lines of vegetation, with a few standalone spots
    % this seems to be independent of drought length -- even a year or two of
    % 0.1 will change our system drastically
    
    % a short-term (<10 years) drought of 0.15 can recover to the mostly
    % vegetated pattern
    
    % a longer term 0.15 drought (say 20 years) gives wider bands of vegetation
    % closer to Klausmeier-type predications
    % r(150:155) = 0.1*r(150:155);
    
    
    % start our record of values, each column is one time step
    Field.biomass_record = zeros(Field.size^2, Initial_Conditions.T+1);
    Field.deep_water_record = zeros(Field.size^2, Initial_Conditions.T+1);
    Field.deep_nitrogen_record = zeros(Field.size^2, Initial_Conditions.T+1);
    % add initial values
    Field.biomass_record(:, 1) = biomass;
    Field.deep_water_record(:, 1) = deep_water;
    Field.deep_nitrogen_record(:, 1) = deep_nitrogen;

    % change water matrix based on slope
    flat_water_matrix = 0.0625*ones(3);
    flat_water_matrix(2, 2) = -0.5;

    if Field.slope_gradient == 0
        Vectors.Water.grass_transport = flat_water_matrix;
    elseif Field.slope_gradient <= 10
        Vectors.Water.grass_transport = flat_water_matrix+(Vectors.Water.grass_transport-flat_water_matrix)*Field.slope_gradient/10;
    else
        Vectors.Water.grass_transport = Vectors.Water.empty_transport;
    end

    % set wind direction
    Vectors.Wind.grass_transport = rot90(Vectors.Wind.grass_transport, Vectors.Wind.direction);
    Vectors.Wind.empty_transport = rot90(Vectors.Wind.empty_transport, Vectors.Wind.direction);
end
