function [ coeffs ] = fit_coefficients( interp, values )

    if strcmpi(interp.type, 'compecon')
        coeffs = funfitxy( interp.fspace, interp.Phi, values ); % what if Phi is not there ?
        
        
    elseif strcmpi( interp.type, 'spinterp')
        spoptions = interp.spoptions;
        
        values_on_grid = values;

        noutput = size(values,2);
        
        surpluses = {};
        
        d = interp.z.d;
        range = interp.z.range;
        s = interp.s;
        
        for i = 1:noutput
            vali = values_on_grid(:,i);
            %bloodydummyfun = @(x) griddatan(s, vali, x); %% VERY SLOW
            bloodydummyfun = @(x) stupidfun(s,vali,x);
            zt = spvals(bloodydummyfun, d, range, spoptions);
            surpluses = [surpluses; zt.vals];
        end
        
        resp = {};
        for i = 1:size(surpluses,1) 
            resp = [ resp ; vertcat( surpluses{i,:} ) ];
        end
        
        coeffs = horzcat(resp{:});
        
    else
        error('Unknown interpolation type');
    end

end

