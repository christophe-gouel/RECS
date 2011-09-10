function [c,x,z,f,exitflag] = recsSolveREEIter(interp,model,s,x,c,options)
% RECSSOLVEREEITER Finds the REE of a model by iteration between equilibrium equations and rational expectations
%
% RECSSOLVEREEITER is called by RECSSOLVEREE. It is not meant to be called directly
% by the user.
%
% See also RECSSOLVEEREE, RECSSOLVEEREEFULL.

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt
  
%% Initialization
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
output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
p      = size(func('h',s(1,:),x(1,:),[],e(1,:),s(1,:),x(1,:),params,output),2);
k      = length(w);               % number of shock values
z      = zeros(n,0);

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
  [c,~,exitflag] = nsoli(@ResidualREE, c(:), reesolveroptions);
  if exitflag==0, exitflag = 1; else exitflag = 0; end

 case 'krylov'
  [c,~,exitflag] = nsoli(@ResidualREE, c(:), reesolveroptions);
  if exitflag==0, exitflag = 1; else exitflag = 0; end

 case 'sa'
  [c,~,exitflag] = SA(@ResidualREE, c(:), reesolveroptions);

 case 'fsolve' % In test - Slow, because it uses numerical derivatives
  if options.display==1
    reesolveroptions = optimset('display','iter-detailed','Diagnostics','on');
  end
  [c,~,exitflag] = fsolve(@ResidualREE, c(:), reesolveroptions);

 case 'kinsol'
  neq = numel(c);
  KINoptions  = KINSetOptions('Verbose',       false,...
                              'LinearSolver',  'GMRES',...
                              'ErrorMessages', false,...
                              'FuncNormTol',   reesolveroptions.atol);
  KINInit(@ResidualREE,neq,KINoptions);
  [status, c] = KINSol(c(:),'LineSearch',ones(neq,1),ones(neq,1));
  KINFree;
  if status==0 || status==1, exitflag = 1; else exitflag = 0; end

end

%%
function [R,FLAG] = ResidualREE(cc)
% RESIDUALREE Calculates the residual of the model with regards to rational expectations

  switch funapprox
    case 'expapprox'
      cc    = reshape(cc,n,p);
      if functional, params{end} = cc; end
      
      z     = funeval(cc,fspace,Phi);
      [x,f] = recsSolveEquilibrium(s,x,z,func,params,cc,e,w,fspace,options);
      
      ind     = (1:n);
      ind     = ind(ones(1,k),:);
      ss      = s(ind,:);
      xx      = x(ind,:);
      output  = struct('F',1,'Js',0,'Jx',0);
      snext   = func('g',ss,xx,[],e(repmat(1:k,1,n),:),[],[],params,output);
      if extrapolate, snextinterp = snext;
      else          
        snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)),...
                          fspace.a(ones(n*k,1),:)); 
      end
      [LB,UB] = func('b',snextinterp,[],[],[],[],[],params);
      
      % xnext calculated by interpolation
      xnext   = min(max(funeval(funfitxy(fspace,Phi,x),fspace,snextinterp),LB),UB);
      if ~useapprox  % xnext calculated by equation solve
        xnext = recsSolveEquilibrium(snext,...
                                     xnext,...
                                     funeval(cc,fspace,snextinterp),...
                                     func,params,cc,e,w,fspace,options);
      end
      
      output              = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
      if nargout(func)<6
        h                 = func('h',ss,xx,[],e(repmat(1:k,1,n),:),snext,xnext,params,output);
      else
        [h,~,~,~,~,hmult] = func('h',ss,xx,[],e(repmat(1:k,1,n),:),snext,xnext,params,output);
        h                 = h.*hmult;
      end
      z     = reshape(w'*reshape(h,k,n*p),n,p);
      
      R     = funfitxy(fspace,Phi,z)-cc;
      
    case 'expfunapprox'
      cc    = reshape(cc,n,p);
      if functional, params{end} = cc; end
      
      [x,f]  = recsSolveEquilibrium(s,x,z,func,params,cc,e,w, fspace,options);
      output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
      R      = funfitxy(fspace,Phi,func('h',[],[],[],[],s,x,params,output))-cc;
      
    case {'resapprox-complete','resapprox-simple'}
      cc    = reshape(cc,n,m);
      if functional, params{end} = cc; end
      
      if useapprox % x calculated by interpolation
        [LB,UB] = func('b',s,[],[],[],[],[],params);
        x       = min(max(funeval(cc,fspace,Phi),LB),UB);
      end % if not previous x is used
      
      if strcmp(funapprox,'resapprox-simple')
        ind    = (1:n);
        ind    = ind(ones(1,k),:);
        ss     = s(ind,:);
        xx     = x(ind,:);
        output = struct('F',1,'Js',0,'Jx',0);
        snext  = func('g',ss,xx,[],e(repmat(1:k,1,n),:),[],[],params,output);
        
        if extrapolate, snextinterp = snext;
        else          
          snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)),...
                            fspace.a(ones(n*k,1),:)); 
        end
        [LB,UB] = func('b',snextinterp,[],[],[],[],[],params);
        xnext   = min(max(funeval(cc,fspace,snextinterp),LB),UB);
        
        output              = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
        if nargout(func)<6
          h                 = func('h',ss,xx,[],e(repmat(1:k,1,n),:),snext,xnext,params,output);
        else
          [h,~,~,~,~,hmult] = func('h',ss,xx,[],e(repmat(1:k,1,n),:),snext,xnext,params,output);
          h                 = h.*hmult;
        end
        z     = reshape(w'*reshape(h,k,n*p),n,p);
      end
      
      [x,f] = recsSolveEquilibrium(s,x,z,func,params,cc,e,w,fspace,options);
      R     = funfitxy(fspace,Phi,x)-cc;
  end
  
  R       = R(:);
  FLAG    = 0;
end

end
