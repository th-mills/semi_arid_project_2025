% function called each time step to create an overall transport matrix
% using the constant matrices 
function sig = generate_sig2(bio, bio_max, unvegetated_transport, vegetated_transport)
    % f is our function for assessing when to use each transport type,
    % depending on cell content
    % we use a row matrix and leave the transport matrices as-is to
    % multiply each column by f(j)
    f = heaviside(bio - 0.1*bio_max)';
    sig = f.*vegetated_transport + (1-f).*unvegetated_transport; 
end
