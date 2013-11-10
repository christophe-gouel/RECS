function [F,J,grid] = recsDeterministicPb(X,fp,gp,hp,s0,xss,p,e,params,nx,grid,nnzJac,izprevT,isprevT,ixT,izT, ...
                                     isT,ixnextT,iznextT,isnextT)
% RECSDETERMINISTICPB Evaluates the equations and Jacobian of the deterministic problem.
%
% RECSDETERMINISTICPB is called by recsSolveDeterministicPb. It is not meant to be
% called directly by the user.
%
% RECSDETERMINISTICPB evaluates the equation defined for each period t by
%
% F = [... f z-h s_{t+1}-g ...];
%
% and its associated Jacobian is defined for each period by
%
%          |       t-1  :      t         :    t+1     |
% -------------------------------------------------------
%          |         s  :   x   z s{t+1} :  x{t+1}    |
% -------------------------------------------------------
% f        |        fs  :  fx  fz      0 :       0    |
% z-h      |       -hs  : -hx  Id   -hsn :    -hxn    |
% s{t+1}-g |       -gs  : -gx   0     Id :       0    |
% ------------------------------------------------------
%          |  Q{t-1,-1} :     Q{t,0}     :  Q{t+1,1}  |
%
% All these per-period Jacobians are arranged to form J, a blocktridiagonal
% matrix (where below T indicates the time horizon):
%
%     ------------------------------------------------------------
%     | Q{0, 0} Q{1, 1}                                          |
%     | Q{0,-1} Q{1, 0} Q{2, 1}                                  |
%     |         Q{1,-1} Q{2, 0} Q{3, 1}                          |
%     |                                 .                        |
% J = |                                   .                      |
%     |                                     .                    |
%     |                              Q{T-2,-1} Q{T-1, 0} Q{T, 1} |
%     |                                        Q{T-1,-1} Q{T, 0} |
%     ------------------------------------------------------------
%
% See also RECSSOLVEDETERMINISTICPB.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
d = size(s0,2);
m = size(xss,2);

n = size(X,1);
M = m+sum(nx);
T = n/(M+p+d);

e = e(ones(T,1),:);

X = reshape(X,M+p+d,T)';

x     = X(:,1:m);
xnext = [x(2:end,:); xss];
w     = X(:,(m+1):(m+nx(1)));
v     = X(:,(m+nx(1)+1):M);
z     = X(:,(M+1):(M+p));
snext = X(:,(M+p+1):(M+p+d));
s     = [s0; snext(1:end-1,:)];

%% Computation of equations and Jacobian
if nargout>=2
  if 1
  %% With Jacobian
  [f,fs,fx,fz]              = fp(s,x,w,v,z,params,ones(4,1));
  [h,hs,hx,~,hsnext,hxnext] = hp(s,x,e,snext,xnext,params,[1 1 1 0 1 1]);
  [g,gs,gx]                 = gp(s,x,e,params,[1 1 1 0]);

  % Main diagonal blocks
  identityz = eye(p);
  identityz = permute(identityz(:,:,ones(1,T)),[3 1 2]);
  identitys = eye(d);
  identitys = permute(identitys(:,:,ones(1,T)),[3 1 2]);
  Jmd = cat(3,...
            cat(2,fx,-hx,-gx),...
            cat(2,fz,identityz,zeros(T,d,p)),...
            cat(2,zeros(T,M,d),-hsnext,identitys));

  % Subdiagonal blocks
  Jsub = cat(3,...
             zeros(T-1,M+p+d,M+p),...
             cat(2,fs(2:end,:,:),-hs(2:end,:,:),-gs(2:end,:,:)));

  % Superdiagonal blocks
  Jsup = cat(3,...
             cat(2,zeros(T-1,M,M),-hxnext(1:end-1,:,:),zeros(T-1,d,M)),...
             zeros(T-1,M+p+d,p+d));

  [J,grid] = blktridiag(permute(Jmd ,[2 3 1]),...
                        permute(Jsub,[2 3 1]),...
                        permute(Jsup,[2 3 1]),...
                        [],grid);
  else
    [f,fs,fx,fz]              = fp(s,x,w,v,z,params,ones(4,1));
    [h,hs,hx,~,hsnext,hxnext] = hp(s,x,e,snext,xnext,params,[1 1 1 0 1 1]);
    [g,gs,gx]                 = gp(s,x,e,params,[1 1 1 0]);

    J                         = spalloc(n,n,nnzJac);
    nd2spdiag = @(Mat) spblkdiag(permute(Mat,[2 3 1]));

    %% df
    J(ixnextT,isprevT) = nd2spdiag(fs(2:end,:,:)); % df/ds
    J(ixT,ixT)         = nd2spdiag(fx);            % df/dx
    J(ixT,izT)         = nd2spdiag(fz);            % df/dz

    %% d(z-h)
    J(iznextT,isprevT)        = nd2spdiag(-hs(2:end,:,:));         % d(z-h)/ds
    J(izT,ixT)                = nd2spdiag(-hx);                    % d(z-h)/dx
    J(sub2ind([n n],izT,izT)) = 1;                                 % d(z-h)/dz
    J(izT,isT)                = nd2spdiag(-hsnext);                % d(z-h)/dsnext
    J(izprevT,ixnextT)        = nd2spdiag(-hxnext(1:(end-1),:,:)); % d(z-h)/dxnext

    %% d(sn-g)
    J(isnextT,isprevT)        = nd2spdiag(-gs(2:end,:,:));         % d(sn-g)/ds
    J(isT,ixT)                = nd2spdiag(-gx);                    % d(sn-g)/dx
    J(sub2ind([n n],isT,isT)) = 1;                                 % d(sn-g)/dsn
  end

else
  %% Without Jacobian
  h = hp(s,x,e,snext,xnext,params,[1 0 0 0 0 0]);
  f = fp(s,x,w,v,z,params,[1 0 0 0]);
  g = gp(s,x,e,params,[1 0 0 0]);

end
F = [f z-h snext-g]';
F = reshape(F,n,1);