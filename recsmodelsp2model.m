function model1 = recsmodelsp2model(model0)

nperiods  = model0.nperiods;

%% Problem's dimensions
dim = model0.dim;

m  = sum(cell2mat(dim(1:nperiods  ,2)))+...
     sum(cell2mat(dim(1:nperiods-1,3)))+...
     sum(cell2mat(dim(2:nperiods  ,1)));
mx = sum(cell2mat(dim(1:nperiods  ,2)));
mz = sum(cell2mat(dim(1:nperiods-1,3)));
ms = sum(cell2mat(dim(2:nperiods  ,1)));

d = dim{1,1};

p = dim{nperiods,3};

E = cell(nperiods,1);
for i=1:nperiods, E{i} = model0.shocks{i}.w'*model0.shocks{i}.e; end

ConstantBounds = model0.bounds;

%%
model1.g = @transition;
model1.h = @expectations;
model1.f = @arbitrage;
model1.b = @bounds;


function [g,gs,gx,ge] = transition(s,x,e,p,o)
  
  n = size(s,1);
  
  xI = x(:,(mx-dim{nperiods,2}):mx);
  sI = x(:,(m-dim{nperiods,1}):end);
    
  [g,gS,gX,ge] = model0.functions(nperiods).g(SI,XI,e,p,o);

  gs = zeros(n,d,d);
  gx = zeros(n,d,m);
  gx(:,:,sum(cell2mat(dim(1:nperiods-1,2))):mx) = gX;
  gx(:,:,(m-dim{nperiods,1}):end)               = gS;
end

function [h,hs,hx,he,hsnext,hxnext] = expectations(s,x,e,snext,xnext,p,o)
  
  sI = x(:,(m-dim{nperiods,1}):end);
  xI = x(:,(mx-dim{nperiods,2}):mx);
  snextI = xnext(:,(m-dim{nperiods,1}):end);
  xnextI = xnext(:,(mx-dim{nperiods,2}):mx);

  [h,hS,hX,he,hsnext,hXNEXT] = model0.functions(nperiods).h(SI,XI,e,SNEXTI,...
                                                    XNEXTI,p,o);
  
  hs = zeros(n,p,d);

  hx = zeros(n,p,m);
  hx(:,:,(mx-dim{nperiods,2}):mx) = hX;
  hx(:,:,(m-dim{nperiods,1}):end) = hS;

  hXNEXT = zeros(n,p,m);
  hXNEXT(:,:,1:dim{1,2}) = hXNEXT;
end

function [f,fs,fx,fz] = arbitrage(s,x,z,p,o)
  
  n = size(s,1);
  
  f  = zeros(n,m);
  fs = zeros(n,m,d);
  fx = zeros(n,m,m);
  fz = zeros(n,m,p);
  
  S = cell(nperiods,1);
  X = cell(nperiods,1);
  Z = cell(nperiods,1);
  S{1} = s;
  Z{nperiod} = z;
  ix = 1;
  for i=1:nperiods
    X{i} = x(:,ix:(ix-1+dim{i,2}))
    ix = ix+dim{i,2};
  end
  for i=1:nperiods-1
    Z{i} = x(:,ix:(ix-1+dim{i,3}))
    ix = ix+dim{i,3};
  end
  for i=2:nperiods
    S{i} = x(:,ix:(ix-1+dim{i,1}))
    ix = ix+dim{i,1};
  end

  %% f
  ix = 1;
  iz = mx+1;
  is = mx+mz+1;
  for i=1:nperiods
    rgf = ix:(ix-1+dim{i,2}); % Range of x
    [f(:,rgf),fS,fX,fZ] = model0.functions(i).f(S{i},X{i},Z{i},p,o);

    if o(1) && (i==1)
      fs(:,rgf,:) = fS;
    end
    if o(2)
      fx(:,rgf,rgf) = fX;
      if i~=nperiods
        fx(:,rgf,iz:(iz-1+dim{i,3})) = fZ;
        iz = iz+dim{i,3};
      end
      if i~=1
        fx(:,rgf,is:(is-1+dim{i,1})) = fS;
        is = is+dim{i,A};
      end
    end
    if o(3) && (i==nperiods)
      fz(:,rgf,:) = fZ;
    end

    ix = ix+dim{i,2};
  end
  
  %% h
  jx = 1;
  is = mx+mz+1;
  for i=1:nperiods-1
    rgf = ix:(ix-1+dim{i,3}); % Range of x
    [f(:,rgf),hS,hX,~,hSNEXT,hXNEXT] = model0.functions(i).h(S{i},X{i},E{i},S{i+1},X{i+1},p,o);

    if o(1) && (i==1)
      fs(:,rgf,:) = hS;
    end
    if o(2)
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
    [f(:,rgf),gS,gX] = model0.functions(i).g(S{i},X{i},E{i},p,o);

    if o(1) && (i==1)
      fs(:,rgf,:) = gS;
    end
    if o(2)
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

function [LB,UB,LB_s,UB_s] = bounds(s,p,o)
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