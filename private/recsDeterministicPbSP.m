function [F,J] = recsDeterministicPbSP(X,functions,s0,xss,e,params,M,P,D,ixT,izT,isT,ixnext,izprev,iznext,nnzJac,dim,istart)
% RECSDETERMINISTICPBSP Evaluates the equations and Jacobian of the deterministic problem with subperiods.
%
% RECSDETERMINISTICPBSP is called by recsSolveDeterministicPbSP. It is not meant
% to be called directly by the user.
%
% RECSDETERMINISTICPBSP evaluates the equation defined for each subperiod i by
%
% F = [... fi zi-hi s_{i+1}-gi ...];
%
% and its associated Jacobian is defined for each subperiod by
%
%           |       i-1   :      i         :    i+1     |
% -------------------------------------------------------
%           |        si   :  xi  zi s{i+1} :  x{i+1}    | 
% -------------------------------------------------------
% fi        |        fs   :  fx  fz      0 :       0    |
% zi-hi     |       -hs   : -hx  Id   -hsn :    -hxn    |
% s{i+1}-gi |       -gs   : -gx   0     Id :       0    |
% ------------------------------------------------------
%           | Q{t,i-1,-1} :    Q{t,i,0}    : Q{t,i+1,1} |
%
% All these per-subperiod Jacobians are arranged to form J, a blocktridiagonal
% matrix (where below i indicates the starting subperiod and T the time
% horizon):
%
%     --------------------------------------------------------------
%     | Q{1,i, 0} Q{1,i+1, 1}                                      |
%     | Q{1,i,-1} Q{1,i+1, 0} Q{1,i+2, 1}                          |
%     |           Q{1,i+1,-1} Q{1,i+2, 0} Q{1,i+3, 1}              |
%     |                                   .                        |
% J = |                                     .                      |
%     |                                       .                    |
%     |                              Q{T,2,-1} Q{T,3, 0} Q{T,4, 1} |
%     |                                        Q{T,3,-1} Q{T,4, 0} |
%     --------------------------------------------------------------
%
% See also RECSSOLVEDETERMINISTICPBSP.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
n        = size(X,1);
T        = length(isT{1});
nperiods = length(e);

inext = @(iperiod) (iperiod+1)*(iperiod<nperiods)+1*(iperiod==nperiods);
X2xzs = @(iX,dimX) reshape(X(iX),dimX,[])';

x       = cellfun(X2xzs,ixT,dim(:,2),'UniformOutput', false);
z       = cellfun(X2xzs,izT,dim(:,3),'UniformOutput', false);
snext   = cellfun(X2xzs,isT,dim(:,1),'UniformOutput', false);
isnextT = isT;
ixnextT = ixT;

isnext = isnextT{inext(istart)}((dim{inext(istart),1}+1):end);

s           = snext;
if istart~=1
  s{istart} = [s0; snext{istart}];
  s{1}      = snext{1}(1:end-1,:);
  isT{1}    = isT{1}(1:end-1,:);
  isprev    = isT{istart};
else
  s{1}      = [s0; snext{1}(1:end-1,:)];
end

xnext           = x;
if istart~=1
  xnext{1}      = [x{1}; xss];
  xnext{istart} = x{istart}(2:end,:);
  ixnextT{istart} = ixnextT{istart}((dim{istart,2}+1):end);
else
  xnext{1}      = [x{1}(2:end,:); xss];
end

f = cell(nperiods,1);
h = cell(nperiods,1);
g = cell(nperiods,1);

%% Computation of equations and Jacobian
if nargout>=2  
  %% With Jacobian
  gs     = cell(nperiods,1);
  gx     = cell(nperiods,1);
  fs     = cell(nperiods,1);
  fx     = cell(nperiods,1);
  fz     = cell(nperiods,1);
  hs     = cell(nperiods,1);
  hx     = cell(nperiods,1);
  hsnext = cell(nperiods,1);
  hxnext = cell(nperiods,1);
  nd2spdiag = @(Mat) spblkdiag(permute(Mat,[2 3 1]));

  %% Calculate basic elements of the equations and Jacobians
  for i=1:nperiods
    % f
    [f{i},fs{i},fx{i},fz{i}] = functions(i).f(s{i},x{i},z{i},params,[1 1 1 1]);

    % z-h
    [h{i},hs{i},hx{i},~,hsnext{i},hxnext{i}] = functions(i).h(s{i},x{i},e{i},...
                                                      snext{inext(i)},...
                                                      xnext{inext(i)},...
                                                      params,[1 1 1 0 1 1]);
    h{i} = z{i}-h{i};
    
    % s-g
    [g{i},gs{i},gx{i}] = functions(i).g(s{i},x{i},e{i},params,[1 1 1 0]);
    g{i} = snext{inext(i)}-g{i};
  end  

  %% Jacobian calculation
  J                         = spalloc(n,n,nnzJac);
  for i=1:nperiods
    %% df
    if i~=istart, J(ixT{i},isT{i})    = nd2spdiag(fs{i});            % df/ds
    else          J(ixnext{i},isprev) = nd2spdiag(fs{i}(2:end,:,:));
    end
    J(ixT{i},ixT{i})                  = nd2spdiag(fx{i});            % df/dx
    J(ixT{i},izT{i})                  = nd2spdiag(fz{i});            % df/dz
    
    %% d(z-h)
    if i~=istart, J(izT{i},isT{i}) = nd2spdiag(-hs{i});                  % d(z-h)/ds
    else          J(iznext,isprev) = nd2spdiag(-hs{i}(2:end,:,:));
    end
    J(izT{i},ixT{i})               = nd2spdiag(-hx{i});                  % d(z-h)/dx
    J(izT{i},isnextT{inext(i)})    = nd2spdiag(-hsnext{i});              % d(z-h)/dsnext
    if i~=nperiods                                                       % d(z-h)/dxnext
      J(izT{i},ixnextT{inext(i)})  = nd2spdiag(-hxnext{i});             
    else
      J(izprev,ixnext{inext(i)}) = nd2spdiag(-hxnext{i}(1:(end-1),:,:));
    end    
    J(sub2ind([n n],izT{i},izT{i})) = 1;                                 % d(z-h)/dz

    %% d(sn-g)
    if i~=istart, J(isnextT{inext(i)},isT{i}) = nd2spdiag(-gs{i});            % d(sn-g)/ds
    else          J(isnext,isprev)            = nd2spdiag(-gs{i}(2:end,:,:));
    end
    J(isnextT{inext(i)},ixT{i}) = nd2spdiag(-gx{i});                          % d(sn-g)/dx
    J(sub2ind([n n],isnextT{i},isnextT{i})) = 1;                              % d(sn-g)/dsn

  end

else
  %% Without Jacobian
  for i=1:nperiods
    f{i} = functions(i).f(s{i},x{i},z{i},params,[1 0 0 0]);
    h{i} = z{i}-functions(i).h(s{i},x{i},e{i},snext{inext(i)},xnext{inext(i)},...
                               params,[1 0 0 0 0 0]);
    g{i} = snext{inext(i)}-functions(i).g(s{i},x{i},e{i},params,[1 0 0 0]);
  
  end

end

%% Reshape output
fstart = cellfun(@(dimX) zeros(0,dimX),dim(:,2),'UniformOutput', false);
hstart = cellfun(@(dimX) zeros(0,dimX),dim(:,1),'UniformOutput', false);
gstart = cellfun(@(dimX) zeros(0,dimX),dim(:,3),'UniformOutput', false);
for i=istart:nperiods
  fstart{i} = f{i}(1,:);
  f{i}      = f{i}(2:end,:);
  hstart{i} = h{i}(1,:);
  h{i}      = h{i}(2:end,:);
  gstart{i} = g{i}(1,:);
  g{i}      = g{i}(2:end,:);
end
Fstart = [fstart hstart gstart]';
Fstart = Fstart(:,istart:end);
F      = [f h g]';
F      = [cat(2,Fstart{:})'; reshape(cat(2,F{:})',(T-1)*(M+D+P),1)];
