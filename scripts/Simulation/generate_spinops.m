function [x y z] = generate_spinops(num_spins, spin_order);
% Feed this a number of spins and the spin order and it gives you the x, y
% and z operators as a matrix of size [2^ns, 2^ns, ns], where x(:, :, n) is
% the Ix operator for the nth spin.
%
% spin_order = 0 gives spin 1/2
% spin_order = 1 gives spin 1
% spin_order = 2 gives spin 3/2

% Set up the Pauli matrices
if spin_order == 0
    sig_x = (1/2)*[0 1; 1 0];
    sig_y = (1/2)*1i*[0 -1; 1 0];
    sig_z = (1/2)*[1 0; 0 -1];
    id = eye(2);
elseif spin_order == 1
    sig_x = (1/sqrt(2))*[0 1 0; 1 0 1; 0 1 0];
    sig_y = (1i/sqrt(2))*[0 -1 0; 1 0 -1; 0 1 0];
    sig_z = [1 0 0; 0 0 0; 0 0 -1];
    id = eye(3);
elseif spin_order == 2
    t = sqrt(3);
    sig_x = (1/2)*[0 t 0 0; t 0 2 0; 0 2 0 t; 0 0 t 0];
    sig_y = (1i/2)*[0 -t 0 0; t 0 -2 0; 0 2 0 -t; 0 0 t 0];
    sig_z = (1/2)*[3 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -3];
    id = eye(4);
end

% Pre-allocate the x, y and z matrices.
ns = num_spins;
expo = spin_order + 2;
x = zeros(expo^ns, expo^ns, ns);
y = zeros(expo^ns, expo^ns, ns);
z = zeros(expo^ns, expo^ns, ns);

for n = 1:ns
    if n == 1
        current_x = sig_x;
        current_y = sig_y;
        current_z = sig_z;
    else
        current_x = id;
        current_y = id;
        current_z = id;
    end
    
    if ns > 1        
        for m = 2:ns
            if m == n
                current_x = kron(current_x, sig_x);
                current_y = kron(current_y, sig_y);
                current_z = kron(current_z, sig_z);
            else
                current_x = kron(current_x, id);
                current_y = kron(current_y, id);
                current_z = kron(current_z, id);
            end
        end
    end
    
    x(:, :, n) = current_x;
    y(:, :, n) = current_y;
    z(:, :, n) = current_z;
end