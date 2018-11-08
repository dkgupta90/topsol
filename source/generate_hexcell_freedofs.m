function [freedofs, fixeddofs] = generate_hexcell_freedofs(ndofx, ndofy, n)

nodofs = [];
fixeddofs = [];
%% generate the tempn for the left side fixeddofs
for i = 1:1:n
    tempk = (1:(n+1-i));
    nodofs = [nodofs ((i-1)*ndofy + tempk)];
    nodofs = [nodofs ((i-1)*ndofy +(ndofy+1))-tempk];
end

%mirroring for the other end
nodofs = [nodofs ((ndofx*ndofy+1)-nodofs)];
freedofs = setdiff(1:ndofy*ndofx, nodofs);

% save the freedofs
save('freedofs_hexcell.dat', 'freedofs', '-ascii');


fixeddofs = [n (ndofy+1-n) ((n-1)*ndofy + 1) (n*ndofy)];
fixeddofs = [fixeddofs ((ndofx*ndofy+1)-fixeddofs)];
end