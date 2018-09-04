element = 'O';
basis = 'aug-cc-pvdz';
spindens = 'beta';

folder = basis;
datafile = strcat(folder, '/', element, '_', spindens);

data = importdata(strcat(datafile,'.dat'));

naos = data.data(1,6);

M = zeros(naos, naos);

i_offset = 1;
for j = 1:naos
    
    j_reduced = rem(j,5); % 1 2 3 4 0 1 2 3 4 0 .... 
    
    if (j_reduced == 0) 
        j_reduced = 5; % Fixes the zero
    end
        
    for i = 1:naos
        
        i_index = i_offset + i;
        
        M(i, j) = data.data(i_index, j_reduced + 1);
       
    end
    
    if (rem(j,5)==0) 
        i_offset = i_offset + naos + 1;
    end
 
end

processed = fopen(strcat(datafile,'.inp'), 'w');
fmt = [repmat('%.16f ', 1, size(M,2)-1), '%.16f\n'];
fprintf(processed, fmt, M);
fclose(processed);

