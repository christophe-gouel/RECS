function [F,J] = pabSimpleEquilibrium(x,s,z,func,params,ggrid,c,e,w,interp)

    % recover x's original size
    [n,d] = size(s);
    m = prod(size(x))/n;

    x     = reshape(x',m,n)';

    K     = size(e,1);
    ind   = (1:n);
    ind   = ind(ones(1,K),:);
    ss    = s(ind,:);   %%%% each observation in ss and xx
    xx    = x(ind,:);   %%%%  is repeated K times
    ee    = e(repmat(1:K,1,n),:);

    if (nargout == 2)  % Analytic derivatives
        [snext,~,gx]           = func('g',ss,xx,[],ee,[],[],params);

        [LB,UB]              = func('b',snext,[],[],[],[],[],params);

        %Xnext = funeval(c,fspace,snext,[zeros(1,d); eye(d)]);  % future controls
        Xnext = eval_function(interp, c, snext, []) ;
        
        xnext                = min(max(Xnext(:,:,1),LB),UB);
        xnextds              = Xnext(:,:,2:end);


        [h,~,hx,hsnext,hxnext] = func('h',ss,xx,[],ee,snext,xnext,params);

        p         = size(h,2);
        z         = reshape(w'*reshape(h,K,n*p),n,p);
        [F,~,fx,fz] = func('f',s,x,z,[],[],[],params);

        F         = reshape(F',n*m,1);

        t1 = hsnext;
        t2 = arraymult(hxnext,xnextds,K*n,p,m,d);  % hx * Cx * gx
        
        Jtmp = hx+arraymult( t1 + t2 ,gx,K*n,p,d,m);   %
        Jtmp = reshape(w'*reshape(Jtmp,K,n*p*m),n,p,m);    %
        
        J    = fx+arraymult(fz,Jtmp,n,m,p,m);
        J    = permute(J,[2 3 1]);
        J = spblkdiag(J,ggrid);

    else % return residual only

        [snext]           = func('g',ss,xx,[],ee,[],[],params);
        [LB,UB]              = func('b',snext,[],[],[],[],[],params);

        %xnext = funeval(c,fspace,snext); % future controls
        xnext = eval_function(interp, c, snext) ;
        xnext                = min(max(xnext,LB),UB);

        h = func('h',ss,xx,[],ee,snext,xnext,params);
        p         = size(h,2);
        z         = reshape(w'*reshape(h,K,n*p),n,p);
        [F] = func('f',s,x,z,[],[],[],params);

        F         = reshape(F',n*m,1);
    end
    
end

