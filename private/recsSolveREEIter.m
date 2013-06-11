function [c,x,z,fval,exitflag] = recsSolveREEIter(interp,model,s,x,c,options)
% RECSSOLVEREEITER Finds the REE of a model by iteration between equilibrium equations and rational expectations
%
% RECSSOLVEREEITER is called by recsSolveREE. It is not meant to be called directly
% by the user.
%
% See also RECSSOLVEEREE, RECSSOLVEREEITERNEWTON, RECSSOLVEEREEFULL.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
explicit           = options.explicit;
extrapolate        = options.extrapolate;
funapprox          = lower(options.funapprox);
functional         = options.functional;
reesolver          = lower(options.reesolver);
reesolveroptions   = catstruct(struct('showiters' , options.display,...
                                      'atol'      , sqrt(eps),...
                                      'lmeth'     , 3,...
                                      'rtol'      , eps),...
                               options.reesolveroptions);
useapprox          = options.useapprox;

b         = model.b;
e         = model.e;
f         = model.f;
g         = model.g;
h         = model.h;
ixforward = model.ixforward;
params    = model.params;
w         = model.w;
[m,p]     = model.dim{2:3};
mf        = sum(ixforward); % Number of forward response variables

fspace = interp.fspace;
Phi    = interp.Phi;

n      = size(x,1);
output = struct('F',1,'Js',0,'Jx',0,'Jz',0,'Jsn',0,'Jxn',0,'hmult',1);
k      = length(w);               % number of shock values
z      = zeros(n,0);

%% Precalculations
if explicit || any(strcmp(funapprox,{'expapprox','resapprox-simple'}))
  [LB,UB] = b(s,params);
  ind     = (1:n);
  ind     = ind(ones(1,k),:);
  ss      = s(ind,:);
  ee      = e(repmat(1:k,1,n),:);
end

if strcmp(funapprox,'resapprox-complete')
  c = c(:,ixforward);
end

%% Solve for the rational expectations equilibrium
switch reesolver
  % Attention: the variables x and z are changed by the nested function 'ResidualREE'
  case 'mixed'
    reesolveroptions.maxit = 10;
    reesolveroptions.atol  = 1E-2;
    reesolveroptions.rtol  = 1E-3;
    c = SA(@ResidualREE, c(:), reesolveroptions);

    reesolveroptions.maxit = 40;
    reesolveroptions.atol  = 1E-7;
    reesolveroptions.rtol  = 1E-25;
    [c,~,exitREE] = nsoli(@ResidualREE, c(:), reesolveroptions);
    if exitREE==0, exitREE = 1; else exitREE = 0; end

  case 'krylov'
    [c,~,exitREE] = nsoli(@ResidualREE, c(:), reesolveroptions);
    if exitREE==0, exitREE = 1; else exitREE = 0; end

  case 'sa'
    [c,~,exitREE] = SA(@ResidualREE, c(:), reesolveroptions);

  case 'kinsol'
    neq = numel(c);
    KINoptions  = KINSetOptions('Verbose',       false,...
                                'LinearSolver',  'GMRES',...
                                'ErrorMessages', false,...
                                'FuncNormTol',   reesolveroptions.atol);
    KINInit(@ResidualREE,neq,KINoptions);
    [status, c] = KINSol(c(:),'LineSearch',ones(neq,1),ones(neq,1));
    KINFree;
    if status==0 || status==1, exitREE = 1; else exitREE = 0; end

end % switch reesolver

if exitREE==1 && exitEQ==1
  exitflag = 1;
else
  exitflag = 0;
end

if strcmp(funapprox,'resapprox-complete')
  c = funfitxy(fspace,Phi,x);
end


%% Nested function
function [R,FLAG] = ResidualREE(cc)
% RESIDUALREE Calculates the residual of the model with regards to rational expectations

  if ~explicit
    switch funapprox
      case 'expapprox'
        %% Expectations approximations
        cc    = reshape(cc,n,p);
        if functional, params{end} = cc; end

        % Calculation of z by interpolation
        z     = funeval(cc,fspace,Phi);

        % Calculation of x
        [x,fval,exitEQ] = recsSolveEquilibrium(s,x,z,b,f,g,h,params,cc,e,w,...
                                               fspace,ixforward,options);

        % Calculation of snext
        xx      = x(ind,:);
        snext   = g(ss,xx,ee,params,output);
        if extrapolate>=1, snextinterp = snext;
        else
          snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)),...
                            fspace.a(ones(n*k,1),:));
        end

        % Calculation of xnext
        [LBnext,UBnext] = b(snext,params);
        if useapprox % xnext calculated by interpolation
          xnext              = zeros(n*k,m);
          xnext(:,ixforward) = min(max(funeval(funfitxy(fspace,Phi,x(:,ixforward)),...
                                               fspace,snextinterp),...
                                       LBnext(:,ixforward)),UBnext(:,ixforward));
        else  % xnext calculated by equation solve
          xnext = min(max(funeval(funfitxy(fspace,Phi,x),fspace,snextinterp),...
                          LBnext),UBnext);
          xnext = recsSolveEquilibrium(snext,...
                                       xnext,...
                                       funeval(cc,fspace,snextinterp),...
                                       b,f,g,h,params,cc,e,w,fspace,...
                                       ixforward,options);
        end

        % Calculation of z
        if nargout(h)<6
          hv                 = h(ss,xx,ee,snext,xnext,params,output);
        else
          [hv,~,~,~,~,hmult] = h(ss,xx,ee,snext,xnext,params,output);
          hv                 = hv.*hmult;
        end
        z     = reshape(w'*reshape(hv,k,n*p),n,p);

        % Prepare output
        R     = funfitxy(fspace,Phi,z)-cc;

      case 'expfunapprox'
        %% Expectations function approximation
        cc    = reshape(cc,n,p);
        if functional, params{end} = cc; end

        [x,fval,exitEQ] = recsSolveEquilibrium(s,x,z,b,f,g,h,params,cc,e,w,...
                                               fspace,ixforward,options);
        output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
        R      = funfitxy(fspace,Phi,h([],[],[],s,x,params,output))-cc;

      case 'resapprox-simple'
        %% Response variables approximation
        cc    = reshape(cc,n,m);
        if functional, params{end} = cc; end

        if useapprox % x calculated by interpolation
          [LB,UB] = b(s,params);
          x       = min(max(funeval(cc,fspace,Phi),LB),UB);
        end % if not previous x is used

        % Calculation of snext
        xx     = x(ind,:);
        snext  = g(ss,xx,ee,params,output);
        
        % Calculation of xnext
        if extrapolate>=1, snextinterp = snext;
        else
          snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)),...
                            fspace.a(ones(n*k,1),:));
        end % extrapolate
        [LBnext,UBnext]    = b(snext,params);
        xnext              = zeros(n*k,m);
        xnext(:,ixforward) = min(max(funeval(cc(:,ixforward),fspace,snextinterp),...
                                     LBnext(:,ixforward)),UBnext(:,ixforward));
        
        % Calculation of z
        if nargout(h)<6
          hv                 = h(ss,xx,ee,snext,xnext,params,output);
        else
          [hv,~,~,~,~,hmult] = h(ss,xx,ee,snext,xnext,params,output);
          hv                 = hv.*hmult;
        end
        z     = reshape(w'*reshape(hv,k,n*p),n,p);

        [x,fval,exitEQ] = recsSolveEquilibrium(s,x,z,b,f,g,h,params,cc,e,w,...
                                               fspace,ixforward,options);
        R            = funfitxy(fspace,Phi,x)-cc;
      
      case 'resapprox-complete'
        %% Response variables approximation
        cc    = reshape(cc,n,mf);
        if functional, params{end} = cc; end

        [x,fval,exitEQ] = recsSolveEquilibrium(s,x,z,b,f,g,h,params,cc,e,w,...
                                               fspace,ixforward,options);
        R            = funfitxy(fspace,Phi,x(:,ixforward))-cc;
    end % switch funapprox
  else
    %% Explicit model
    cc = reshape(cc,n,m);

    % Calculation of x by interpolation
    if useapprox, x = min(max(funeval(cc,fspace,Phi),LB),UB); end

    % Calculation of snext
    xx     = x(ind,:);
    snext  = g(ss,xx,ee,params,output);

    % Calculation of xnext
    if extrapolate>=1, snextinterp = snext;
    else
      snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)),...
                        fspace.a(ones(n*k,1),:));
    end % extrapolate
    [LBnext,UBnext]    = b(snext,params);
    xnext              = zeros(n*k,m);
    xnext(:,ixforward) = min(max(funeval(cc(:,ixforward),fspace,snextinterp),...
                                 LBnext(:,ixforward)),UBnext(:,ixforward));

    % Calculation of z
    if nargout(h)<6
      hv                 = h(ss,xx,ee,snext,xnext,params,output);
    else
      [hv,~,~,~,~,hmult] = h(ss,xx,ee,snext,xnext,params,output);
      hv                 = hv.*hmult;
    end
    z      = reshape(w'*reshape(hv,k,n*p),n,p);

    % Calculation of x by explicit formula
    x      = min(max(f(s,[],z,params,output),LB),UB);

    % Prepare output
    R      = funfitxy(fspace,Phi,x)-cc;
    exitEQ = 1;
    fval      = zeros(size(x));

  end % if ~explicit
  R        = R(:);
  FLAG     = 0;
end

end
