function final_graphs(Initial_Conditions, Field, grass, rain, maps, tick_size, time_diff)
%FINAL_GRAPHS Summary of this function goes here
%   Detailed explanation goes here

    % output final biomass as a larger plot for visibility
    figure(Name = 'Final Biomass', NumberTitle = 'off');
    biomass_out = Field.biomass_record(:, Initial_Conditions.T+1);
    biomass_out = reshape(biomass_out, Field.size, Field.size);
    biomass_out = biomass_out(2:Field.size, 1:Field.size);
    imagesc(biomass_out)
    colormap(maps.grassmap)
    colorbar
    % clim([0 grass.b_max])
    axis square
    axis ij
    xticks(0:tick_size:Field.size)
    yticks(0:tick_size:Field.size)
    title("Final Biomass")

    % output mean biomass against time, along with rainfall for comparison
    mean_biomass = mean(Field.biomass_record)';
    figure(Name = 'Mean Biomass', NumberTitle = 'off');
    plot(mean_biomass, 'LineWidth', 1)
    title("Mean Biomass against Time")
    xlim([0 Initial_Conditions.T])
    xticks(0:time_diff:Initial_Conditions.T)
    xlabel("Time (years)")
    ylim([0 grass.b_max])
    yticks(0:30:grass.b_max)
    ylabel("Mean Biomass (g/m2)")
    % add horizontal line to see convergence
    yline(mean_biomass(Initial_Conditions.T+1,1), "--")
    
    figure(Name = 'Rainfall', NumberTitle = 'off');
    plot(rain, 'LineWidth', 1)
    title("Rainfall")
    xlim([0 Initial_Conditions.T])
    xticks(0:time_diff:Initial_Conditions.T)
    xlabel("Time (years)")
    ylim([0 450])
    ylabel("Yearly Rainfall (g/m2)")
    
    % output mean soil resources against time, on one pair of axes
    mean_water = mean(Field.deep_water_record)'; 
    mean_nitrogen = mean(Field.deep_nitrogen_record)';
    figure(Name = 'Soil Resources', NumberTitle = 'off');
    title("Soil Resources against Time")
    xlabel("Time (years)")
    xlim([0 Initial_Conditions.T])
    
    yyaxis left
    plot(mean_water, 'LineWidth', 1)
    ylabel("Mean Soil Water (g/m2)")
    % ylim([0 water_saturation])
    % yline(mean_water(Initial_Conditions.T+1, 1), "--")
    
    yyaxis right
    plot(mean_nitrogen, 'LineWidth', 1)
    ylabel("Mean Soil Nitrogen (g/m2)")
    % ylim([0 nitrogen_saturation])
    % yline(mean_nitrogen(Initial_Conditions.T+1, 1), ":")

    % code for generating and playing movie of biomass
    % doesn't work particularly well
    frames = Initial_Conditions.T/5;
    tstep = Initial_Conditions.T/frames;
    M(frames+1) = struct('cdata', [], 'colormap', maps.grassmap);
    fig = figure(Name = 'Movie', NumberTitle = 'off');
    axis square
    axis ij
    xlim([1 Field.size-1])
    ylim([1 Field.size-1])
    xticks(0:tick_size:Field.size)
    yticks(0:tick_size:Field.size)
    colormap(maps.grassmap)
    colorbar
    clim([0 200])
    ax = gca;
    ax.NextPlot = 'replaceChildren';
    title 'Biomass over Time'
    for j=1:frames+1
        movie_out = Field.biomass_record(:, tstep*(j-1)+1);
        movie_out = reshape(movie_out, Field.size, Field.size);
        movie_out = movie_out(2:Field.size, 1:Field.size);
        imagesc(ax, movie_out)
        xlabel('t = ' + string(tstep*(j-1)))
        M(j) = getframe(fig);
    end
    framerate = 5;
    movie(fig, M, 4, framerate);