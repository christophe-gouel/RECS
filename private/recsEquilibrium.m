function [F,Jx,Jc] = recsEquilibrium(x,s,z,b,f,g,h,params,gridJx,c,e,w,fspace,funapprox,extrapolate,ixforward,ArrayProblem)
% RECSEQUILIBRIUM evaluates the equilibrium equations and Jacobian
%
% RECSEQUILIBRIUM is called by recsSolveREEFull and recsSolveEquilibrium. It is not
% meant to be called directly by the user.
%
% See also RECSSOLVEREEFULL, RECSSOLVEEQUILIBRIUM.

% Copyright (C) 2011-2016 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
n     = size(s,1);
if ~ArrayProblem, x = reshape(x,[],n)'; end
m     = size(x,2);
mnext = length(ixforward); % Number of next-period response variables
mf    = sum(ixforward); % Number of forward response variables

%% Evaluate the equilibrium equations and Jacobian
switch funapprox
  case 'expapprox'
    if nargout==2
      %% With Jacobian
      [F,~,Jx] = f(s,x,z,params,[1 0 1 0]);
    else
      %% Without Jacobian
      F        = f(s,x,z,params);
    end
  case {'expfunapprox','resapprox'}
    k     = size(e,1);
    ind   = (1:n);
    ind   = ind(ones(1,k),:);
    ss    = s(ind,:);
    xx    = x(ind,:);
    ee    = e(repmat(1:k,1,n),:);

    if nargout>=2
      %% With Jacobians
      [snext,~,gx]        = g(ss,xx,ee,params,[1 0 1 0]);
      d = size(snext,2);
      if extrapolate>=1, snextinterp = snext;
      else
          snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)),fspace.a(ones(n*k,1),:));
      end
      Bsnext = funbasx(fspace,snextinterp,[zeros(1,d); eye(d)],'expanded');
      % It seems to be faster with 'expanded', but may use more memory. To confirm
      % later.
%       Bsnext = funbasx(fspace,snextinterp,[zeros(1,d); eye(d)]);

      switch funapprox
        case 'expfunapprox'
          H                 = funeval(c,fspace,Bsnext,[zeros(1,d); eye(d)]);
          hv = H(:,:,1);
          hs = H(:,:,2:end);
        case 'resapprox'
          [LBnext,UBnext]      = b(snext,params);
          Xnext                = funeval(c,fspace,Bsnext,...
                                         [zeros(1,d); eye(d)]);
          xnext                = zeros(n*k,mnext);
          xnext(:,ixforward)   = min(max(Xnext(:,:,1),LBnext(:,ixforward)),...
                                     UBnext(:,ixforward));
          xfnextds             = Xnext(:,:,2:end);

          [hv,~,hx,~,hsnext,hxnext] = h(ss,xx,ee,snext,xnext,params,[1 0 1 0 1 1]);
      end
      p           = size(hv,2);
      z           = reshape(w'*reshape(hv,k,n*p),n,p);
      [F,~,fx,fz] = f(s,x,z,params,[1 0 1 1]);

      switch funapprox
        case 'expfunapprox'
          Jxtmp = arraymult(hs,gx,k*n,p,d,m);
        case 'resapprox'
          Jxtmp = hx+...
                  arraymult(hsnext+...
                            arraymult(hxnext(:,:,ixforward),xfnextds,k*n,p,mf,d),...
                            gx,k*n,p,d,m);
      end
      Jxtmp = reshape(w'*reshape(Jxtmp,k,n*p*m),n,p,m);
      Jx    = fx+arraymult(fz,Jxtmp,n,m,p,m);

      if nargout==3
        %% With Jacobian with respect to c
        if ~strcmp(Bsnext.format,'expanded')
          Bsnext = funbconv(Bsnext,zeros(1,d));
        end
        Bsnext = mat2cell(Bsnext.vals{1}',n,k*ones(n,1))';
        switch funapprox % The product with fz is vectorized for 'expfunapprox' and
                         % executed in a loop for resapprox', because it
                         % seems to be the fastest ways to do it.
          case 'expfunapprox'
            if issparse(Bsnext{1})
              Jc = cellfun(@(X) kron((X*w)',speye(p)),Bsnext,'UniformOutput',false);
              Jc = cat(1,Jc{:});
              fz = spblkdiag(permute(fz,[2 3 1]));
              Jc = fz*Jc;        
            else
              Jc = cellfun(@(X) kron((X*w)',  eye(p)),Bsnext,'UniformOutput',false);
              Jc = permute(reshape(cat(1,Jc{:}),[p n numel(c)]),[2 1 3]);
              Jc = arraymult(fz,Jc,n,m,p,numel(c));
              Jc = reshape(permute(Jc,[2 1 3]),[n*m numel(c)]);
            end
            
          case 'resapprox'
            [~,gridJc] = spblkdiag(zeros(p,mf,k),[],0);
            if issparse(Bsnext{1})
              kw     = kron(w',speye(p));
              hxnext = num2cell(reshape(hxnext(:,:,ixforward),[k n p mf]),[1 3 4])';
              Jctmp  = cellfun(...
                  @(X,Y) kw*spblkdiag(permute(X,[3 4 2 1]),gridJc)*kron(Y',speye(mf)),...
                  hxnext,Bsnext,'UniformOutput',false);
              fz = spblkdiag(permute(fz,[2 3 1]));
              Jc = fz*cat(1,Jctmp{:});
            else
              kw     = kron(w',eye(p));
              hxnext = num2cell(reshape(hxnext(:,:,ixforward),[k n p mf]),[1 3 4])';
              Jctmp  = cellfun(...
                  @(X,Y) kw*full(spblkdiag(permute(X,[3 4 2 1]),gridJc))*kron(Y',eye(mf)),...
                  hxnext,Bsnext,'UniformOutput',false);
              Jc   = zeros(n*m,numel(c));
              for i=1:n
                Jc((i-1)*m+1:i*m,:) = permute(fz(i,:,:),[2 3 1])*Jctmp{i};
              end
            end
        end % funapprox
      end
    else
      %% Without Jacobian
      snext   = g(ss,xx,ee,params);

      switch funapprox
        case 'expfunapprox'
          if extrapolate>=1, snextinterp = snext;
          else
            snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)), ...
                              fspace.a(ones(n*k,1),:));
          end
          hv                  = funeval(c,fspace,snextinterp);

        case 'resapprox'
          [LBnext,UBnext] = b(snext,params);
          if extrapolate>=1, snextinterp = snext;
          else
            snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)), ...
                              fspace.a(ones(n*k,1),:));
          end
          xnext              = zeros(n*k,mnext);
          xnext(:,ixforward) = min(max(funeval(c,fspace,snextinterp),...
                                       LBnext(:,ixforward)),UBnext(:,ixforward));
          hv                 = h(ss,xx,ee,snext,xnext,params);
      end
      p       = size(hv,2);
      z       = reshape(w'*reshape(hv,k,n*p),n,p);
      F       = f(s,x,z,params);
    end
end

% Reshape output
if ~ArrayProblem
  F = reshape(F',n*m,1);
  if nargout>=2
    Jx = permute(Jx,[2 3 1]);
    if n>1, Jx = spblkdiag(Jx,gridJx); end
  end
end
