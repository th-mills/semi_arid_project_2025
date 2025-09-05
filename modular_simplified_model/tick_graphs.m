function tick_graphs(Initial_Conditions, Field, maps, tick_size, t)
%TICK_GRAPHS Summary of this function goes here
%   Detailed explanation goes here

    % get current output row vectors and convert to N x N matrices
    biomass_out = Field.biomass_record(:,t);
    biomass_out = reshape(biomass_out, Field.size, Field.size);
    % cutting out the top line, which is overly vegetated due to boundary
    % conditions
    biomass_out = biomass_out(2:Field.size, 1:Field.size);

    deep_water_out = Field.deep_water_record(:,t);
    deep_water_out = reshape(deep_water_out, Field.size, Field.size);

    deep_nitrogen_out = Field.deep_nitrogen_record(:,t);
    deep_nitrogen_out = reshape(deep_nitrogen_out, Field.size, Field.size);

    % create a figure showing biomass and soil resources at each step
    figure(Name = 't = ' + string(Initial_Conditions.start_year+t-1), NumberTitle = 'off');
    tiledlayout(1,3)

    ax1 = nexttile;
    imagesc(biomass_out)
    colormap(ax1, maps.grassmap)
    colorbar
    % clim([0 Grass.b_max])
    axis square
    axis ij
    xticks(0:tick_size:Field.size)
    yticks(0:tick_size:Field.size)
    title(ax1, "Biomass at t = " + string(Initial_Conditions.start_year + t-1))

    ax2 = nexttile;
    imagesc(deep_water_out)
    colormap(ax2, maps.watermap)
    colorbar
    % clim([0 water_saturation])
    axis square
    axis ij
    xticks(0:tick_size:Field.size)
    yticks(0:tick_size:Field.size)
    title(ax2, "Soil Water at t = " + string(Initial_Conditions.start_year + t-1))

    ax3 = nexttile;
    imagesc(deep_nitrogen_out)
    colormap(ax3, maps.nitrogenmap)
    colorbar
    % clim([0 nitrogen_saturation])
    axis square
    axis ij
    xticks(0:tick_size:Field.size)
    yticks(0:tick_size:Field.size)
    title(ax3, "Soil Nitrogen at t = " + string(Initial_Conditions.start_year + t-1))
