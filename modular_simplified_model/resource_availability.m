% function for calculating availability
function av = resource_availability(bio, res, N, bio_max, periodic)
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

    if periodic
        % fudge attempt 1: average out top and bottom rows' availability 
        average = (av(1,:) + av(N,:))./2;
        av(1,:) = average;
        av(N,:) = average;
    end
    
    av = reshape(av, N^2, 1);
end
