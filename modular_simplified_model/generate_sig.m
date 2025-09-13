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





% ---------- OLD VERSION OF ABOVE FUNCTION ------------------------------ %

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
%}
