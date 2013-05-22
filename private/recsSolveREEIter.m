function [c,x,z,f,exitflag] = recsSolveREEIter(interp,model,s,x,c,options)
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

e      = model.e;
func   = model.func;
params = model.params;
w      = model.w;

fspace = interp.fspace;
Phi    = interp.Phi;

[n,m]  = size(x);
output = struct('F',1,'Js',0,'Jx',0,'Jz',0,'Jsn',0,'Jxn',0,'hmult',1);
p      = size(func('h',s(1,:),x(1,:),[],e(1,:),s(1,:),x(1,:),params,output),2);
k      = length(w);               % number of shock values
z      = zeros(n,0);

%% Precalculations
if explicit || any(strcmp(funapprox,{'expapprox','resapprox-simple'}))
  [LB,UB] = func('b',s,[],[],[],[],[],params);
  ind     = (1:n);
  ind     = ind(ones(1,k),:);
  ss      = s(ind,:);
  ee      = e(repmat(1:k,1,n),:);
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
        [x,f,exitEQ] = recsSolveEquilibrium(s,x,z,func,params,cc,e,w,fspace,options);

        % Calculation of snext
        xx      = x(ind,:);
        snext   = func('g',ss,xx,[],ee,[],[],params,output);
        if extrapolate>=1, snextinterp = snext;
        else
          snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)),...
                            fspace.a(ones(n*k,1),:));
        end

        % Calculation of xnext
        [LBnext,UBnext] = func('b',snext,[],[],[],[],[],params);
        % xnext calculated by interpolation
        xnext   = min(max(funeval(funfitxy(fspace,Phi,x),fspace,snextinterp),LBnext),UBnext);
        if ~useapprox  % xnext calculated by equation solve
          xnext = recsSolveEquilibrium(snext,...
                                       xnext,...
                                       funeval(cc,fspace,snextinterp),...
                                       func,params,cc,e,w,fspace,options);
        end

        % Calculation of z
        if nargout(func)<6
          h                 = func('h',ss,xx,[],ee,snext,xnext,params,output);
        else
          [h,~,~,~,~,hmult] = func('h',ss,xx,[],ee,snext,xnext,params,output);
          h                 = h.*hmult;
        end
        z     = reshape(w'*reshape(h,k,n*p),n,p);

        % Prepare output
        R     = funfitxy(fspace,Phi,z)-cc;

      case 'expfunapprox'
        %% Expectations function approximation
        cc    = reshape(cc,n,p);
        if functional, params{end} = cc; end

        [x,f,exitEQ] = recsSolveEquilibrium(s,x,z,func,params,cc,e,w, fspace,options);
        output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
        R      = funfitxy(fspace,Phi,func('h',[],[],[],[],s,x,params,output))-cc;

      case {'resapprox-complete','resapprox-simple'}
        %% Response variables approximation
        cc    = reshape(cc,n,m);
        if functional, params{end} = cc; end

        if useapprox % x calculated by interpolation
          [LB,UB] = func('b',s,[],[],[],[],[],params);
          x       = min(max(funeval(cc,fspace,Phi),LB),UB);
        end % if not previous x is used

        if strcmp(funapprox,'resapprox-simple')
          % Calculation of snext
          xx     = x(ind,:);
          snext  = func('g',ss,xx,[],ee,[],[],params,output);

          % Calculation of xnext
          if extrapolate>=1, snextinterp = snext;
          else
            snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)),...
                              fspace.a(ones(n*k,1),:));
          end % extrapolate
          [LBnext,UBnext] = func('b',snext,[],[],[],[],[],params);
          xnext   = min(max(funeval(cc,fspace,snextinterp),LBnext),UBnext);

          % Calculation of z
          if nargout(func)<6
            h                 = func('h',ss,xx,[],ee,snext,xnext,params,output);
          else
            [h,~,~,~,~,hmult] = func('h',ss,xx,[],ee,snext,xnext,params,output);
            h                 = h.*hmult;
          end
          z     = reshape(w'*reshape(h,k,n*p),n,p);
        end % if strcmp(funapprox,'resapprox-simple')

        [x,f,exitEQ] = recsSolveEquilibrium(s,x,z,func,params,cc,e,w,fspace,options);
        R            = funfitxy(fspace,Phi,x)-cc;
    end % switch funapprox
  else
    %% Explicit model
    cc = reshape(cc,n,m);

    % Calculation of x by interpolation
    if useapprox, x = min(max(funeval(cc,fspace,Phi),LB),UB); end

    % Calculation of snext
    xx     = x(ind,:);
    snext  = func('g',ss,xx,[],ee,[],[],params,output);

    % Calculation of xnext
    if extrapolate>=1, snextinterp = snext;
    else
      snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)),...
                        fspace.a(ones(n*k,1),:));
    end % extrapolate
    [LBnext,UBnext] = func('b',snext,[],[],[],[],[],params);
    xnext   = min(max(funeval(cc,fspace,snextinterp),LBnext),UBnext);

    % Calculation of z
    if nargout(func)<6
      h                 = func('h',ss,xx,[],ee,snext,xnext,params,output);
    else
      [h,~,~,~,~,hmult] = func('h',ss,xx,[],ee,snext,xnext,params,output);
      h                 = h.*hmult;
    end
    z      = reshape(w'*reshape(h,k,n*p),n,p);

    % Calculation of x by explicit formula
    x      = min(max(func('f',s,[],z,[],[],[],params,output),LB),UB);

    % Prepare output
    R      = funfitxy(fspace,Phi,x)-cc;
    exitEQ = 1;
    f      = zeros(size(x));

  end % if ~explicit
  R        = R(:);
  FLAG     = 0;
end

end
