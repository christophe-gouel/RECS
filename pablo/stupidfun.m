function [res] = stupidfun(x, y, z)

    n = size(x,1);
    d = size(x,2);
    
    q = size(y,2);
    
    nz = size(z,1);
    
    res = zeros(nz, q);
    
    for i = 1:nz
       ind = 0;
       for j = 1:n
           if z(i,:) == x(j,:);
               res(i,:) = y(j,:);
               ind = j;
               break;
           end
       end
       if ind == 0;
           res(i,:) = nan;
       end
    end

end