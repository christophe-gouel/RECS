function [ allres ] = eval_function( interp, coeffs, points, derivs )

    if strcmpi(interp.type, 'compecon')
        fspace = interp.fspace;
        if nargin == 3;
            allres = funeval(coeffs, fspace, points);
        elseif nargin == 4;
            d = interp.fspace.d;
            allres = funeval(coeffs, fspace, points,[zeros(1,d); eye(d)]);
        end
    elseif strcmpi(interp.type, 'spinterp')
        
        
        n = size(coeffs,1);
        d = size(coeffs,2);
        % painful but needed : split coffs into arrays
        
        template = interp.z.vals;
        
        % vals will contain the hierarchical surpluses to be evaluated
        vals = cell( d, size(template,2) );
        t0 = 1;
        for j = 1: size(vals,2)
            nl =  size(template{1,j},1)  ;% number of coefficients to copy        
            for i = 1:d
                vals{i,j} = coeffs(t0:(t0+nl-1),i);
            end
            t0 = t0 + nl;            
        end
        
        z = interp.z;
        z.vals = vals; % we shouldn't be changing it every time...
        noutput = size(coeffs,2);
        
        %s = points;
        %x = s;

        if nargin == 3
            with_derivatives = 0;
        elseif nargin == 4
            with_derivatives = 1;
        end

        n = size(points,1);
           
        if with_derivatives
    
            %allres = zeros(n, noutput);
            allres = zeros(n, noutput,1+d);
    
            for i = 1:noutput
                z.selectOutput = i;
                [res,dres] = spinterp(z, points);
                dres = horzcat(dres{:})';
                allres(:,i,1) = res;
                allres(:,i,2:end) = dres;
            end
    
        else
    
            allres = zeros(n, noutput);
    
            for i = 1:noutput
                z.selectOutput = i;
                [res,dres] = spinterp(z, points);
                dres = horzcat(dres{:})';
                allres(:,i) = res;
            end
    
        end
        
    elseif strcmpi(interp.type, 'spderivcc')
        
    else
        error('Unknown interpolation type')
    end


end

