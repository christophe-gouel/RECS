function [ interp ] = create_sp_grid(order, l, gridtype, smin, smax  )

    range = [smin', smax'];

    spoptions = spset;
    spoptions.GridType = gridtype;
    
    spoptions.KeepGrid = 'on';
    spoptions.Vectorized = 'off';
    
    spoptions.SparseIndices = 'off';
    
    spoptions.MinDepth = l;
    spoptions.MaxDepth = l;
    spoptions.FunctionArgType = 'vector';
    
    dumbdummyfun = @(x) x(:,1)*0; % used to compute the grid
    
    zzz = spvals( dumbdummyfun , order, range, spoptions);
    
    % do not obfuscate user's mind
    
    zzz = rmfield(zzz, {'estRelError','estAbsError','fevalRange','minGridVal','maxGridVal','fevalTime','surplusCompTime'});
    
    interp = struct;
    interp.z = zzz;
    interp.s = cell2mat( zzz.grid' );
    interp.type = 'spinterp';
    
    spoptions.KeepGrid = 'off';
    %spoptions.Vectorized = 'on';
    interp.spoptions = spoptions; % keep options
        
end

