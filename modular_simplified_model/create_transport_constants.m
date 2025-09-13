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





% ---------- OLD VERSION OF THE ABOVE FUNCTION -------------------------- %

%{
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
