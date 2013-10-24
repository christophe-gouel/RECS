function model1 = recsmodelsp2model(model0)

nperiods  = model0.nperiods;

%% Problem's dimensions
dim = model0.dim;

m  = sum(cell2mat(dim(1:nperiods  ,2)))+...
     sum(cell2mat(dim(1:nperiods-1,3)))+...
     sum(cell2mat(dim(2:nperiods  ,1)));
mx = sum(cell2mat(dim(1:nperiods  ,2)));
mz = sum(cell2mat(dim(1:nperiods-1,3)));
ms = sum(cell2mat(dim(2:nperiods  ,1))); %#ok<NASGU>

d = dim{1,1};

p = dim{nperiods,3};

E = cell(nperiods,1);
for ip=1:nperiods, E{ip} = model0.shocks{ip}.w'*model0.shocks{ip}.e; end
q = size(E{nperiods},2);

ConstantBounds = model0.bounds;

%%
model1.functions.b  = @bounds;
model1.functions.ee = @(s,x,z,params) NaN;
model1.functions.f  = @arbitrage;
model1.functions.g  = @transition;
model1.functions.h  = @expectations;

model1.functions.bp = @bounds;
model1.functions.fp = @(s,x,w,v,z,params,o) arbitrage(s,x,z,params,o);
model1.functions.gp = @transition;
model1.functions.hp = @expectations;

model1.dim = {d,m,p,q};

model1.params = model0.params;

model1.shocks = model0.shocks{nperiods};

model1.infos.eq_type     = 'mcp';
model1.infos.ixforward   = true(1,m);
model1.infos.ixvarbounds = false(m,2);
model1.infos.model_type  = 'fgh2';
model1.infos.nxvarbounds = zeros(1,2,'int16');

function [g,gs,gx,ge] = transition(s,x,e,params,o)
  
  if nargin <= 4
    o = [1; zeros(3,1)];
    for i = 2:nargout
      o(i) = 1;
    end
  end

  n = size(s,1);
  
  xI = x(:,(mx+1-dim{nperiods,2}):mx);
  sI = x(:,(m+1-dim{nperiods,1}):end);
    
  [g,gS,gX,ge] = model0.functions(nperiods).g(sI,xI,e,params,o);

  gs = zeros(n,d,d);
  gx = zeros(n,d,m);
  if o(3)
    gx(:,:,sum(cell2mat(dim(1:nperiods-1,2))):mx) = gX;
    gx(:,:,(m+1-dim{nperiods,1}):end)               = gS;
  end
end

function [h,hs,hx,he,hsnext,hxnext] = expectations(s,x,e,~,xnext,params,o)
  
  if nargin <= 6
    o = [1; zeros(5,1)];
    for i = 2:nargout
      o(i) = 1;
    end
  end

  n = size(s,1);
  hs = zeros(n,p,d);
  hx = zeros(n,p,m);
  hxnext = zeros(n,p,m);
  
  sI = x(:,(m+1-dim{nperiods,1}):end);
  xI = x(:,(mx+1-dim{nperiods,2}):mx);
  snextI = xnext(:,(m+1-dim{nperiods,1}):end);
  xnextI = xnext(:,(mx+1-dim{nperiods,2}):mx);

  [h,hS,hX,he,hsnext,hXNEXT] = model0.functions(nperiods).h(sI,xI,e,snextI,...
                                                    xnextI,params,o);

  if o(3)
    hx(:,:,(mx+1-dim{nperiods,2}):mx) = hX;
    hx(:,:,(m+1-dim{nperiods,1}):end) = hS;
  end
  if o(6), hxnext(:,:,1:dim{1,2}) = hXNEXT; end
end

function [f,fs,fx,fz] = arbitrage(s,x,z,params,o)
  
  if nargin <= 4
    o = [1; zeros(3,1)];
    for i = 2:nargout
      o(i) = 1;
    end
  end

  n = size(s,1);
  
  f  = zeros(n,m);
  fs = zeros(n,m,d);
  fx = zeros(n,m,m);
  fz = zeros(n,m,p);
  
  S = cell(nperiods,1);
  X = cell(nperiods,1);
  Z = cell(nperiods,1);
  S{1} = s;
  Z{nperiods} = z;
  ix = 1;
  for i=1:nperiods
    X{i} = x(:,ix:(ix-1+dim{i,2}));
    ix = ix+dim{i,2};
  end
  for i=1:nperiods-1
    Z{i} = x(:,ix:(ix-1+dim{i,3}));
    ix = ix+dim{i,3};
  end
  for i=2:nperiods
    S{i} = x(:,ix:(ix-1+dim{i,1}));
    ix = ix+dim{i,1};
  end

  %% f
  ix = 1;
  iz = mx+1;
  is = mx+mz+1;
  for i=1:nperiods
    rgf = ix:(ix-1+dim{i,2}); % Range of x
    [f(:,rgf),fS,fX,fZ] = model0.functions(i).f(S{i},X{i},Z{i},params,o);

    if o(2) && (i==1)
      fs(:,rgf,:) = fS;
    end
    if o(3)
      fx(:,rgf,rgf) = fX;
      if i~=nperiods
        fx(:,rgf,iz:(iz-1+dim{i,3})) = fZ;
        iz = iz+dim{i,3};
      end
      if i~=1
        fx(:,rgf,is:(is-1+dim{i,1})) = fS;
        is = is+dim{i,1};
      end
    end
    if o(4) && (i==nperiods)
      fz(:,rgf,:) = fZ;
    end

    ix = ix+dim{i,2};
  end
  
  %% h
  jx = 1;
  is = mx+mz+1;
  for i=1:nperiods-1
    rgf = ix:(ix-1+dim{i,3}); % Range of x
    [f(:,rgf),hS,hX,~,hSNEXT,hXNEXT] = model0.functions(i).h(S{i},X{i},E{i},S{i+1},X{i+1},params,...
                                                      [1 1 1 0 1 1]);

    if o(2) && (i==1)
      fs(:,rgf,:) = hS;
    end
    if o(3)
      fx(:,rgf,jx:(jx-1+dim{i,2})) = fx(:,rgf,jx:(jx-1+dim{i,2}))+hX;
      fx(:,rgf,jx+dim{i,2}:(jx-1+dim{i,2}+dim{i+1,2})) = fx(:,rgf,jx+dim{i,2}:(jx-1+dim{i,2}+dim{i+1,2}))+hXNEXT;

      if i~=1
        fx(:,rgf,is:(is-1+dim{i,1})) = fx(:,rgf,is:(is-1+dim{i,1}))+hS;
      end
      fx(:,rgf,(is+dim{i,1}):(is-1+dim{i,1}+dim{i+1,1})) = fx(:,rgf,(is+dim{i,1}):(is-1+dim{i,1}+dim{i+1,1}))+hSNEXT;
      is = is+dim{i,1};
    end
    
    ix = ix+dim{i,3};
    jx = jx+dim{i,2};
  end
  
  %% g
  jx = 1;
  is = mx+mz+1;
  for i=1:nperiods-1
    rgf = ix:(ix-1+dim{i+1,1}); % Range of f
    [f(:,rgf),gS,gX] = model0.functions(i).g(S{i},X{i},E{i},params,o);

    if o(2) && (i==1)
      fs(:,rgf,:) = gS;
    end
    if o(3)
      fx(:,rgf,jx:(jx-1+dim{i,2})) = gX;
      if i~=1
        fx(:,rgf,is:(is-1+dim{i,1})) = gS;
      end
      is = is+dim{i,1};
    end
    
    ix = ix+dim{i,1};
    jx = jx+dim{i,2};
  end
  
end

function [LB,UB,LB_s,UB_s] = bounds(s,~,~)
% This transformation is only valid for constant bounds.
  
  n = size(s,1);
  
  LB = -inf(n,m);
  UB =  inf(n,m);
  LB_s = zeros(n,m,d);
  UB_s = zeros(n,m,d);
  
  ix = 1;
  for i=1:nperiods
    rgx = ix:(ix-1+dim{i,2}); % Range of x
    LB(:,rgx) = repmat(ConstantBounds{i,1},n,1);
    UB(:,rgx) = repmat(ConstantBounds{i,2},n,1);
    ix = ix+dim{i,2};
  end
end

end