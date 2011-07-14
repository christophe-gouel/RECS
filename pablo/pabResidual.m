function [R,FLAG] = pabResidual(cc,n,m,p,s,model,interp,options)


%fspace = interp.fspace;
%Phi = interp.Phi;

func = str2func(model.func);
params = model.params;

e = model.e;
w = model.w;


[LB,UB] = func('b',s,[],[],[],[],[],params);

xx = eval_function( interp, cc, s );

x0 = reshape(xx,n,m);

x0 = min(max(x0,LB),UB);   % impose bounds if they are not already satisfied

%     c = funfitxy(fspace,Phi,x0);

% c = cc;  % future controls


% now we look for x, such that recsEquilibrium(s,x,z,func,params,grid,c,e,w,fspace,method) = 0

[~,ggrid] = spblkdiag(zeros(m,m,n),[],0); % computed once for all to save computational time

x       = reshape(x0',[n*m 1]);
LB      = reshape(LB',[n*m 1]);
UB      = reshape(UB',[n*m 1]);

fobj = @(X) pabSimpleEquilibrium( X, s,[],func,params,ggrid,cc,e,w,interp);

eqsolveroptions = options.eqsolveroptions;

[x,f,exitflag] = lmmcp(fobj,x,LB,UB,eqsolveroptions);

if exitflag~=1, disp('No convergence'); end


x    = reshape(x,m,n)';

c = fit_coefficients(interp, x);


%    R     = funfitxy(fspace,Phi,x)-c;
R = c - cc;

z     = [];

% R       = R(:);
FLAG    = 0;

end
