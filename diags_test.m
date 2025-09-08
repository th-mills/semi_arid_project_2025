clear   
% leave at 2800 or below! slows down unbelievably at 2900 for some reason
N = 6;


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
        transport(i, i+down+left) = local_matrix(1,3);
        transport(i, i+down+right) = local_matrix(1,1);
    end
end

function transport = create_transport_constants2(N, local_matrix)
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

    % add left, middle, and right columns' movement to the relevant
    % diagonals
    transport = transport + kron(central_diagonal, middle_block);
    transport = transport + kron(above_diagonal, left_block);
    transport = transport + kron(below_diagonal, right_block);
    
    % add left-right periodicity with a couple of extra blocks
    transport((N^2-N+1):N^2, 1:N) = left_block;
    transport(1:N, (N^2-N+1):N^2) = right_block;
end

% for regular diffusion, we can create D (equivalent to the blocks in the
% above function)
id = speye(N, N);
D = spdiags([1 -2 1], -1:1, N, N);
tensor_lacplacian = kron(D, id) + kron(id, D);
spy(kron(D, id))

% for the downward diagonal motion of the Stewart et al. transport
% matrices, something similar-looking can be done using a matrix of ones
% just below the diagonal
id_down = spdiags(1, -1, N, N);
D2 = spdiags([1 0 1], -1:1, N ,N);
tensor_downslope = kron(D2, id_down) - 2*kron(id, id);

% checking that we get similar output from the function, as with the
% Kronecker product approach
laplacian_matrix = [0 1 0; 1 -4 1; 0 1 0];
function_laplacian = create_transport_constants2(N, laplacian_matrix);
downslope_matrix = [0 0 0; 0 -2 0; 1 0 1];
function_downslope = create_transport_constants2(N, downslope_matrix);

model_movement = [1 4 7; 2 5 8; 3 6 9];
% model_movement = [0 1 0; 1 -5 1; 0.5 1 0.5];
function_total_movement = create_transport_constants2(N, model_movement);

figure;
spy(tensor_lacplacian)
title('Global Laplacian')
for k=0:N-1
    xline(k*N+0.5, ':')
    yline(k*N+0.5, ':')
end

figure;
spy(function_laplacian)
title('Global Laplacian with Periodicity')
for k=0:N-1
    xline(k*N+0.5, ':')
    yline(k*N+0.5, ':')
end

figure;
spy(tensor_downslope)
title('Downslope by Tensors')
for k=0:N-1
    xline(k*N+0.5, ':')
    yline(k*N+0.5, ':')
end

figure;
spy(function_downslope)
title('Downslope with Periodicity')
for k=0:N-1
    xline(k*N+0.5, ':')
    yline(k*N+0.5, ':')
end

figure;
spy(function_total_movement)
title('Movement for Stewart et al. Grass')
for k=0:N-1
    xline(k*N+0.5, ':')
    yline(k*N+0.5, ':')
end
