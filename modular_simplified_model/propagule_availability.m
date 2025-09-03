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