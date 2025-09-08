% INPUTS

% we work on an N x N grid of cells
N = 100;

% model runs for T time steps
T = 600;

% maximum biomass for grass in a cell 
b_max = 319;

% constants for soil saturation approach
water_saturation = 100;
nitrogen_saturation = 5;

% value for which type of rain data is used:
% 0 for average of 243, 1 for historical, 2 for repeated historical, 
% leaving 3 for potential stochastic 
% rainfall is generated in 'Initialise Model' section below
rain_type = 2;

% value for type of initial condition used:
% 0 for uniform random distribution, 1 for stripes, 2 for spots, 3 for one
% central square of biomass, 4 for for biomass around the rim
initial_biomass_type = 0;

% likewise for nitrogen
d(1, 1:T) = 1.5;

% paper has:
% w_m = 0.7; n_m = 0.125; w_e = 3.5; n_e = 0.62
% given inputs to Stewart et al. model are:
% w_m = 3.5; n_m = 0.0273; w_e = 17.5; n_e = 0.0273
% using the values from the model

% add constants for plant behaviour
water_maintenance = 3.5;
nitrogen_maintenance = 0.125;
water_efficiency = 17.5;
nitrogen_efficiency = 0.62;
% k is a total non-resource-related mortality rate, i.e. from diseases,
% grazing, physical damage
k = 0.1;
% f is failure rate for propagules to establish themselves
f = 0.05;
% max_growth is maximum growth rate for grass in a year
% paper lists 0.125, while model takes 1.125
% this is likely a difference in notation, since Creosote bush similarly
% has 0.09 in paper and 1.09 in model input
% however the model code appears to use 1.125*biomass as a limit for
% propagules, so this one is used
max_growth = 1.125;

% transport matrices
grass_transport = [0 0.1 0; 0.1 -0.5 0.1; 0.05 0.1 0.05];
empty_transport = [0 0 0; 0 -1 0; 0 1 0];


% DEFINE FUNCTIONS

% function for calculating availability
function av = resource_availability(bio, res, N, bio_max)
    % reshaping into a matrix gives a simpler function, since we are moving
    % downslope and downslope elements would be separated in the vector
    b = reshape(bio, N, N);

    % we can take either a single number for uniform resources, or a vector
    % for cell-specific redistribution of propagules
    % (use a heaviside here to ensure that 1-bio/bio_max doesn't give a
    % small negative due to floating point errors)
    % av = res.*spdiags(heaviside(1-bio./bio_max), 0, N^2, N^2)*(1-bio./bio_max);
    av = res.*heaviside(1-bio./bio_max).*(1-bio./bio_max);
    % change this to an N x N matrix for simpler iteration
    av = reshape(av, N, N);

    % the proportion of the available water moving downslope is given by
    % 0.1 for vegetated cells and 1 for unvegetated cells
    sigs = 1 - 0.9.*heaviside(b - 0.1*bio_max);
    for i=2:N
        % water is carried downslope, including the previously updated
        % values above it
        av(i,:) = av(i,:) + sigs(i-1,:).*av(i-1,:);
    % output availability as a vector
    end
    av = reshape(av, N^2, 1);
end

function av = propagule_availability(bio, prop, N, bio_max)
    % function for propagule redistribution, very similar to resource
    % availability
    b = reshape(bio, N, N);
    % propagules are available proportional to biomass, so highly vegetated
    % spread propagules rather than collecting them
    av = prop.*(bio./bio_max);
    av = reshape(av, N, N);
    sigs = 1 - 0.9.*heaviside(b - 0.1*bio_max);
    for i=2:N
        av(i,:) = av(i,:) + sigs(i-1,:).*av(i-1,:);
    end
    av = reshape(av, N^2, 1);
end

% old functions for local transport, replaced by faster ones below
%{
function sig = generate_sig(bio, N, bio_max)
    % create empty sparse matrix with space for our (up to) seven nonzero 
    % values per row
    sig = spalloc(N^2, N^2, 7*N^2);
    % f is N^2-dim vector  
    f = heaviside(bio - 0.1*bio_max);
    for i=1:N^2
        up = -1;
        left = -N;
        right = +N;
        down = +1;
        % if at top or bottom, periodic conditions apply
        if mod(i, N) == 1
            up = N-1;
        elseif mod(i, N) == 0 
            down = 1-N;
        end
        % if at left or right border, periodic conditions apply
        if i <= N
            left = N*(N-1); 
        elseif i > N*(N-1)
            right = -N*(N-1);
        end
        % add the relevant edge values to the matrix
        sig(i, i) = -1 + 0.5*f(i);
        sig(i, i+up) = 1 - 0.9*f(i+up);
        sig(i, i+down) = 0.1*f(i+down);
        sig(i, i+left) = 0.1*f(i+left);
        sig(i, i+right) = 0.1*f(i+right);
        sig(i, i+up+left) = 0.05*f(i+up+left);
        sig(i, i+up+right) = 0.05*f(i+up+right);
        % water doesn't have any motion diagonally upwards, so the
        % down-left and down-right cells do not contribute to i 

        % also note this is for maximum gradient relevant to our parameters
        % any higher (above 10 in Stewart et al.) just has everything flow
        % downhill without the local transport near biomass
        % if gradient is less, then transport is linearly interpolated with
        % uniform spread: 
        % [0.0625 0.0625 0.0625
        %  0.0625  0.5   0.0625
        %  0.0625 0.0625 0.0625]
    end
end

% matrix to create initial N^2 x N^2 local transport matrices, which are
% then applied to cells by componentwise heaviside 
function transport = create_transport_constants(N, local_matrix)
    % take a 3x3 local transport matrix as our input
    transport = spalloc(N^2, N^2, 7*N^2);
    for i=1:N^2
        up = -1;
        left = -N;
        right = +N;
        down = +1;
        % if at top or bottom, periodic conditions apply
        if mod(i, N) == 1
            up = N-1;
        elseif mod(i, N) == 0 
            down = 1-N;
        end
        % if at left or right border, periodic conditions apply
        if i <= N
            left = N*(N-1); 
        elseif i > N*(N-1)
            right = -N*(N-1);
        end
        % add the relevant edge values to the matrix
        % i + ... is the origin for the flow, while i is the destination
        % so entries of matrix are 'flipped' in both directions
        transport(i, i) = local_matrix(2,2);
        transport(i, i+up) = local_matrix(3,2);
        transport(i, i+down) = local_matrix(1,2);
        transport(i, i+left) = local_matrix(2,3);
        transport(i, i+right) = local_matrix(2,1);
        transport(i, i+up+left) = local_matrix(3,3);
        transport(i, i+up+right) = local_matrix(3,1);
    end
end

%}

% function to generate large constant matrices for local transport before
% the main simulation loop; works for any 3x3 local matrix
function transport = create_transport_constants(N, local_matrix)
    % create matrices with ones along diagonal, and diagonal lines just
    % above and below true diagonal
    central_diagonal = speye(N, N);
    above_diagonal = spdiags(1, 1, N, N);
    below_diagonal = spdiags(1, -1, N, N);

    % create block matrices for movement given by each column of
    % local_matrix
    left_block = spdiags(flip(local_matrix(:, 1))', -1:1, N, N);
    middle_block = spdiags(flip(local_matrix(:, 2))', -1:1, N, N);
    right_block = spdiags(flip(local_matrix(:, 3))', -1:1, N, N);
    
    % add vertical periodicity to each of these
    left_block(N, 1) = local_matrix(1, 1);
    left_block(1, N) = local_matrix(3, 1);
    middle_block(N, 1) = local_matrix(1, 2);
    middle_block(1, N) = local_matrix(3, 2);
    right_block(N, 1) = local_matrix(1, 3);
    right_block(1, N) = local_matrix(3, 3);
    
    % create empty sparse matrix with space for full local transport
    % matrices (i.e. 9 connections for each of our N^2 cells in the region)
    transport = spalloc(N^2, N^2, 9*N^2);

    transport = transport + kron(central_diagonal, middle_block);
    transport = transport + kron(above_diagonal, left_block);
    transport = transport + kron(below_diagonal, right_block);
    
    % add left-right periodicity with a couple of extra blocks
    transport((N^2-N+1):N^2, 1:N) = left_block;
    transport(1:N, (N^2-N+1):N^2) = right_block;
end

% function called each time step to create an overall transport matrix
% using the constant matrices 
function sig = generate_sig(bio, bio_max, unvegetated_transport, vegetated_transport)
    % f is our function for assessing when to use each transport type,
    % depending on cell content
    % we use a row matrix and leave the transport matrices as-is to
    % multiply each column by f(j)
    f = heaviside(bio - 0.1*bio_max)';
    sig = f.*vegetated_transport + (1-f).*unvegetated_transport; 
end

% INITIALISE MODEL

% generate initial biomass distribution
switch initial_biomass_type
    case 0
        % uniform random
        biomass = b_max*0.3*rand(N^2,1);

    case 1
        % stripes 
        biomass = 0.35*b_max*rand(N, N);
        for i=1:N
            biomass(i,:) = 0.5*(1+sin(16*pi*i/N))*biomass(i,:);
        end
        biomass = reshape(biomass, N^2,1);

    case 2
        % spots
        biomass = 0.6*b_max*rand(N, N);
        for i=1:N
            biomass(i,:) = 0.5*(1+sin(10*pi*i/N))*biomass(i,:);
        end
        for j=1:N
            biomass(:,j) = 0.5*(1+sin(10*pi*j/N))*biomass(:,j);
        end
        biomass = reshape(biomass, N^2, 1);

    case 3
        % central spot
        biomass = zeros(N);
        biomass(0.4*N:0.6*N, 0.4*N:0.6*N) = 0.5*b_max*rand(0.2*N+1);
        biomass = reshape(biomass, N^2, 1);

    case 4
        % biomass around sides
        biomass = 0.3*b_max*rand(N);
        biomass(0.2*N:0.8*N, 0.2*N:0.8*N) = 0;
        biomass = reshape(biomass, N^2, 1);
end

% generate initial soil resource distribution
% used 25 and 0.5 from initial values given in Stewart et al. 
deep_water = zeros(N^2, 1);
deep_water(1:N^2, 1) = 25;
deep_nitrogen = zeros(N^2, 1);
deep_nitrogen(1:N^2, 1) = 0.5;

% rainfall uses an average of 243 mm/year by default
r = zeros(T, 1);
r(1:T, 1) = 243;
start_year = 0;
switch rain_type
    % historical rainfall
    case 1
    % get the rain data as a matrix, with years in first column and
    % rainfall in second
    raindata = readmatrix("raindat.dat");
    % change our start year to that of the historical data
    start_year = raindata(1, 1);
    if T > size(raindata, 1)
        % if our simulation goes beyond the length of the record, change
        % only those up to the time-length of the record
        r(1:size(raindata, 1)) = raindata(:, 2);
    else
        r(1:T) = raindata(1:T, 2);
    end
    
    % repeated historical rainfall
    case 2
        raindata = readmatrix("raindat.dat");
        % change our start year to that of the historical data
        start_year = raindata(1, 1);
        len = size(raindata, 1);
        % get quotient and remainder of T / raindata length
        remainder = mod(T, len);
        quotient = (T-remainder)/len;
        % add repeated historical record
        for n=0:(quotient-1)
            r(n*len+1:(n+1)*len) = raindata(:, 2);
        end
        % add final section which is less than the full raindata length
        r(quotient*len+1:quotient*len+remainder) = raindata(1:remainder, 2); 
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
biomass_record = zeros(N^2, T+1);
deep_water_record = zeros(N^2, T+1);
deep_nitrogen_record = zeros(N^2, T+1);
% add initial values
biomass_record(:, 1) = biomass;
deep_water_record(:, 1) = deep_water;
deep_nitrogen_record(:, 1) = deep_nitrogen;

empty_transport_constant = create_transport_constants(N, empty_transport);
grass_transport_constant = create_transport_constants(N, grass_transport);

% MAIN
for t = 1:T

    % produce availability vectors for water and nitrogen
    water_availability = resource_availability(biomass, r(t), N, b_max);
    nitrogen_availability = resource_availability(biomass, d(t), N, b_max);

    % generate matrix for local transport of resources and propagules
    sigma = generate_sig(biomass, b_max, empty_transport_constant, grass_transport_constant);

    % find final surface resource distributions using local transport and
    % availability
    surface_water = (r(t) + sigma*water_availability);
    surface_nitrogen = d(t) + sigma*nitrogen_availability;

    % calculate WR and NR for use in biomass calculations
    water_res = (surface_water + deep_water - water_maintenance.*biomass)./water_efficiency;
    nitrogen_res = (surface_nitrogen + deep_nitrogen + nitrogen_maintenance.*biomass)./nitrogen_efficiency;

    % we construct a 3xN^2 matrix where each row is the resource remaining or maximum
    % growth possible over the cells, by concatenating row vectors
    % min returns the minimum for each column, i.e. the 'I' value for each cell
    % as a row vector, which we then change back to a column vector
    intermediate = min([water_res'; nitrogen_res'; max_growth.*biomass'])';

    % any positive I values mean propagules are produced by a cell, otherwise
    % we have zero biomass in a cell or insufficient resources and none are produced
    propagules = heaviside(intermediate).*intermediate;

    % calculate availability from this
    prop_availability = propagule_availability(biomass, propagules, N, b_max);
        
    % calculate insufficiency terms for the biomass updating equation
    % (this is not in the written form, but makes code more readable)
    water_insufficiency = (1-k)*(water_efficiency/water_maintenance).*(1-heaviside(intermediate)).*heaviside(nitrogen_res-water_res).*intermediate;
    nitrogen_insufficiency = (1-k)*(nitrogen_efficiency/nitrogen_maintenance).*(1-heaviside(intermediate)).*heaviside(water_res-nitrogen_res).*intermediate;
    
    % output some values for the central cell, used to investigate deep
    % soil accumulating too much
    % disp('At time t = ' + string(t+1) + ': B = ' + string(biomass(1250)) + '; W = ' + string(surface_water(1250)) + '; DW = ' + string(deep_water(1250)) + '; I = ' + string(intermediate(1250)))

    % update our soil resource stores 
    % if our cell is nearly empty, we remove resources to prevent
    % accumulation with:
    % (1-0.95.*heaviside(0.1*b_max - biomass)).*
    deep_water = water_efficiency.*(1-0.95.*heaviside(0.1*b_max - biomass)).*(water_res - intermediate);
    deep_nitrogen = nitrogen_efficiency.*(1-0.95.*heaviside(0.1*b_max - biomass)).*(nitrogen_res - intermediate);
    
    % manually prevent too much resource from accumulating 
    % brute-force solution to resources accumulating too much in cells
    deep_water = water_saturation.*heaviside(deep_water - water_saturation) + heaviside(water_saturation - deep_water).*deep_water;
    deep_nitrogen = nitrogen_saturation.*heaviside(deep_nitrogen - nitrogen_saturation) + heaviside(nitrogen_saturation - deep_nitrogen).*deep_nitrogen;

    % update biomass
    biomass = (1-k)*biomass + (1-f)*propagules + (1-f)*(sigma*prop_availability) + water_insufficiency + nitrogen_insufficiency;

    % if any biomass is above the maximum, reduce it to that
    % (note that floating point errors are protected against inside the
    % availability function, so we can set directly to maximum, rather than
    % just below it)
    biomass = b_max.*heaviside(biomass-b_max) + heaviside(b_max-biomass).*biomass;

    % ensure that biomass is non-negative
    % (seen a couple instances of a single cell spiralling to large negative
    % values. may be caused by memory problem or a float error from biomass being near zero)
    biomass = heaviside(biomass).*biomass;

    % add values to the record array as a new column
    biomass_record(:, t+1) = biomass;
    deep_water_record(:, t+1) = deep_water;
    deep_nitrogen_record(:, t+1) = deep_nitrogen;

end

% OUTPUTS
grassmap = uint8([240 200 60; 230 200 60; 220 200 60; 210 200 60; 200 200 60; 180 200 60; 160 200 60; 140 200 60; 120 200 60]);
watermap = uint8([60 200 240; 50 180 220; 40 160 200; 30 140 180; 20 120 180; 10 100 160; 0 80 140]);
nitrogenmap = uint8([50 30 20; 80 40 18; 110 50 16; 140 60 14; 170 70 12; 200 80 10; 230 90 10]);

% calculate an appropriate tick size for axes
tick_size = (N - mod(N, 5))/5;
if tick_size == 0
    tick_size = 1;
end

% find appropriate separation for ten timed outputs
time_diff = (T - mod(T, 10))/10;
if time_diff == 0
    time_diff = 1;
end

for t=1:time_diff:T+1
    % get current output row vectors and convert to N x N matrices
    biomass_out = biomass_record(:,t);
    biomass_out = reshape(biomass_out, N, N);
    % cutting out the top line, which is overly vegetated due to boundary
    % conditions
    biomass_out = biomass_out(2:N, 1:N);

    deep_water_out = deep_water_record(:,t);
    deep_water_out = reshape(deep_water_out, N, N);

    deep_nitrogen_out = deep_nitrogen_record(:,t);
    deep_nitrogen_out = reshape(deep_nitrogen_out, N, N);

    % create a figure showing biomass and soil resources at each step
    figure(Name = 't = ' + string(start_year+t-1), NumberTitle = 'off');
    tiledlayout(1,3)

    ax1 = nexttile;
    imagesc(biomass_out)
    colormap(ax1, grassmap)
    colorbar
    % clim([0 200])
    axis square
    axis ij
    xticks(0:tick_size:N)
    yticks(0:tick_size:N)
    title(ax1, "Biomass at t = " + string(start_year + t-1))

    ax2 = nexttile;
    imagesc(deep_water_out)
    colormap(ax2, watermap)
    colorbar
    % clim([0 water_saturation])
    axis square
    axis ij
    xticks(0:tick_size:N)
    yticks(0:tick_size:N)
    title(ax2, "Soil Water at t = " + string(start_year + t-1))

    ax3 = nexttile;
    imagesc(deep_nitrogen_out)
    colormap(ax3, nitrogenmap)
    colorbar
    % clim([0 nitrogen_saturation])
    axis square
    axis ij
    xticks(0:tick_size:N)
    yticks(0:tick_size:N)
    title(ax3, "Soil Nitrogen at t = " + string(start_year + t-1))
end

% output final biomass as a larger plot for visibility
figure(Name = 'Final Biomass', NumberTitle = 'off');
biomass_out = biomass_record(:, T+1);
biomass_out = reshape(biomass_out, N, N);
biomass_out = biomass_out(2:N, 1:N);
imagesc(biomass_out)
colormap(grassmap)
colorbar
% clim([0 b_max])
axis square
axis ij
xticks(0:tick_size:N)
yticks(0:tick_size:N)
title("Final Biomass")

% output mean biomass against time, along with rainfall for comparison
mean_biomass = mean(biomass_record)';
figure(Name = 'Mean Biomass', NumberTitle = 'off');
plot(mean_biomass, 'LineWidth', 1)
title("Mean Biomass against Time")
xlim([0 T])
xticks(0:time_diff:T)
xlabel("Time (years)")
ylim([0 b_max])
yticks(0:30:b_max)
ylabel("Mean Biomass (g/m2)")
% add horizontal line to see convergence
yline(mean_biomass(T+1,1), "--")

figure(Name = 'Rainfall', NumberTitle = 'off');
plot(r, 'LineWidth', 1)
title("Rainfall")
xlim([0 T])
xticks(0:time_diff:T)
xlabel("Time (years)")
ylim([0 450])
ylabel("Yearly Rainfall (g/m2)")

% output mean soil resources against time, on one pair of axes
mean_water = mean(deep_water_record)'; 
mean_nitrogen = mean(deep_nitrogen_record)';
figure(Name = 'Soil Resources', NumberTitle = 'off');
title("Soil Resources against Time")
xlabel("Time (years)")
xlim([0 T])

yyaxis left
plot(mean_water, 'LineWidth', 1)
ylabel("Mean Soil Water (g/m2)")
% ylim([0 water_saturation])
% yline(mean_water(T+1, 1), "--")

yyaxis right
plot(mean_nitrogen, 'LineWidth', 1)
ylabel("Mean Soil Nitrogen (g/m2)")
% ylim([0 nitrogen_saturation])
% yline(mean_nitrogen(T+1, 1), ":")


% code for generating and playing movie of biomass
% doesn't work particularly well

%{
frames = T/5;
tstep = T/frames;
M(frames+1) = struct('cdata', [], 'colormap', grassmap);
fig = figure(Name = 'Movie', NumberTitle = 'off');
axis square
axis ij
xlim([1 N-1])
ylim([1 N-1])
xticks(0:tick_size:N)
yticks(0:tick_size:N)
colormap(grassmap)
colorbar
clim([0 200])
ax = gca;
ax.NextPlot = 'replaceChildren';
title 'Biomass over Time'
for j=1:frames+1
    movie_out = biomass_record(:, tstep*(j-1)+1);
    movie_out = reshape(movie_out, N, N);
    movie_out = movie_out(2:N, 1:N);
    imagesc(ax, movie_out)
    xlabel('t = ' + string(tstep*(j-1)))
    M(j) = getframe(fig);
end
framerate = 5;
movie(fig, M, 4, framerate);
%}