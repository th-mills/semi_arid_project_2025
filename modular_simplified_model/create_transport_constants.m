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
