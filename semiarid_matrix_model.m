% INPUTS

% we work on an N x N grid of cells
N = 50;

% model runs for T time steps
T = 500;

% maximum biomass for grass in a cell 
b_max = 319;

% constants for soil saturation approach
water_saturation = 100;
nitrogen_saturation = 15;

% value for which type of rain data is used:
% 0 for average of 243, 1 for historical, leaving 2 for potential
% stochastic 
% rainfall is generated in 'Initialise Model' section below
rain_type = 1;

% likewise for nitrogen
d(1, 1:T) = 1.5;

% add constants for plant behaviour
water_maintenance = 0.7;
nitrogen_maintenance = 0.125;
water_efficiency = 3.5;
nitrogen_efficiency = 0.62;
% k is a total non-resource-related mortality rate, i.e. from diseases,
% grazing, physical damage
k = 0.1;
% f is failure rate for propagules to establish themselves
f = 0.05;
% max_growth is maximum growth rate for grass in a year
max_growth = 0.125;

% DEFINE FUNCTIONS

% function for calculating availability
function av = availability(bio, res, N, bio_max)
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
    end
end

% INITIALISE MODEL

% generate initial biomass distribution
biomass = b_max*0.3*rand(N^2,1);
% generate initial soil resource distribution
% used 25 and 0.5 from initial values given in Stewart et al. 
deep_water(1:N^2, 1) = 25;
deep_nitrogen(1:N^2, 1) = 0.5;

% rainfall uses an average of 243 mm/year by default
r = zeros(T, 1);
r(1:T, 1) = 243;
start_year = 0;
if rain_type == 1
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
end

% manually add droughts (fun to play with!)

% 0.1 gives sparse lines of vegetation, with a few standalone spots
% this seems to be independent of drought length -- even a year or two of
% 0.1 will change our system drastically

% a short-term (<10 years) drought of 0.15 can recover to the mostly
% vegetated pattern

% a longer term 0.15 drought (say 20 years) gives wider bands of vegetation
% closer to Klausmeier-type predications
r(350:351) = 0.1*243;


% start our record of values, each column is one time step
biomass_record = zeros(N^2, T+1);
deep_water_record = zeros(N^2, T+1);
deep_nitrogen_record = zeros(N^2, T+1);
% add initial values
biomass_record(:, 1) = biomass;
deep_water_record(:, 1) = deep_water;
deep_nitrogen_record(:, 1) = deep_nitrogen;

% MAIN
for t = 1:T

    % produce availability vectors for water and nitrogen
    water_availability = availability(biomass, r(t), N, b_max);
    nitrogen_availability = availability(biomass, d(t), N, b_max);

    % generate matrix for local transport of resources and propagules
    sigma = generate_sig(biomass, N, b_max);

    % find final surface resource distributions using local transport and
    % availability
    surface_water = r(t) + sigma*water_availability;
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
    propagule_availability = availability(biomass, propagules, N, b_max);
        
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
    biomass = (1-k)*biomass + (1-f)*propagules + (1-f)*(sigma*propagule_availability) + water_insufficiency + nitrogen_insufficiency;

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

    deep_water_out = deep_water_record(:,t);
    deep_water_out = reshape(deep_water_out, N, N);

    deep_nitrogen_out = deep_nitrogen_record(:,t);
    deep_nitrogen_out = reshape(deep_nitrogen_out, N, N);

    % create a figure showing biomass and soil resources at each step
    figure;
    tiledlayout(1,3)

    ax1 = nexttile;
    imagesc(biomass_out)
    colormap(ax1, grassmap)
    colorbar
    clim([0 b_max])
    axis square
    axis ij
    xticks(0:tick_size:N)
    yticks(0:tick_size:N)
    title(ax1, "Biomass at t = " + string(start_year + t-1))

    ax2 = nexttile;
    imagesc(deep_water_out)
    colormap(ax2, watermap)
    colorbar
    clim([0 water_saturation])
    axis square
    axis ij
    xticks(0:tick_size:N)
    yticks(0:tick_size:N)
    title(ax2, "Soil Water at t = " + string(start_year + t-1))

    ax3 = nexttile;
    imagesc(deep_nitrogen_out)
    colormap(ax3, nitrogenmap)
    colorbar
    clim([0 nitrogen_saturation])
    axis square
    axis ij
    xticks(0:tick_size:N)
    yticks(0:tick_size:N)
    title(ax3, "Soil Nitrogen at t = " + string(start_year + t-1))
end

% output final biomass as a larger plot for visibility
figure;
imagesc(biomass_out)
colormap(grassmap)
colorbar
clim([0 b_max])
axis square
axis ij
xticks(0:tick_size:N)
yticks(0:tick_size:N)
title("Final Biomass")

% output mean biomass against time, along with rainfall for comparison
mean_biomass = mean(biomass_record)';
figure;
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

figure;
disp(r)
plot(r, 'LineWidth', 1)
title("Rainfall")
xlim([0 T])
xticks(0:time_diff:T)
xlabel("Time (years)")
ylim([0 400])
ylabel("Yearly Rainfall (g/m2)")

% output mean soil resources against time, on one pair of axes
mean_water = mean(deep_water_record)'; 
mean_nitrogen = mean(deep_nitrogen_record)';
figure;
title("Soil Resources against Time")
xlabel("Time (years)")
xlim([0 T])

yyaxis left
plot(mean_water, 'LineWidth', 1)
ylabel("Mean Soil Water (g/m2)")
ylim([0 water_saturation])
% yline(mean_water(T+1, 1), "--")

yyaxis right
plot(mean_nitrogen, 'LineWidth', 1)
ylabel("Mean Soil Nitrogen (g/m2)")
ylim([0 nitrogen_saturation])
% yline(mean_nitrogen(T+1, 1), ":")