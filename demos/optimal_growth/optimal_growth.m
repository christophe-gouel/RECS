%close all

clear interp model options
clear mex

 warning('off', 'catstruct:DuplicatesFound')

%dbstop if error


% COMPUTE SHOCK DISTRIBUTION
model = struct;
mm = [0.0]';
sigma1 = 0.01;


K = [10];                                           % number of shocks
%K = 3
[model.e,model.w] = qnwnorm(K,mm,sigma1);                    % shocks and proabilities


% PACK MODEL STRUCTURE
model.func = 'optimal_growth_recs';                          % model function file
model.params =  feval(model.func,'params');     % other parameters


% DEFINE APPROXIMATION SPACE
interp = struct;
order     = [ 5 50  ];                             % degree of approximation
smin  = [    (1-sigma1*1.5)  0.5 ];     % minimum state
smax  = [    (1+sigma1*1.5)  5.0 ];     % maximum state
interp.fspace = fundefn('spli',order,smin,smax);       % function space
snodes = funnode(interp.fspace);                   % state collocaton nodes
s=gridmake(snodes);
n             = prod(order);
interp.Phi    = funbasx(interp.fspace);



%% initial guess

dsol = feval(model.func,'solution'); % deterministic solution
[LB,UB] = feval(model.func,'b',s); % bounds

s_ss = repmat(dsol.s_ss,n,1);
x_ss = repmat(dsol.x_ss,n,1);
xinit = x_ss + (s - s_ss) * dsol.Z';

xinit = max(min(xinit,UB-0.001),LB+0.001);
%xinit = max(min(xinit,UB),LB);

%%%%%%%%%%%%%%%%%%%%%%%%%

% compute initialization 

x = xinit;

K = length(model.w);
[n,m] = size(x);

%output = struct('F',true)
p     = size(feval(model.func,'h',s(1,:),x(1,:),[],model.e(1,:),s(1,:),x(1,:),model.params),2);
ind   = (1:n);
ind   = ind(ones(1,K),:);
ss    = s(ind,:);
xx    = x(ind,:);
snext = feval(model.func,'g',ss,xx,[],model.e(repmat(1:K,1,n),:),[],[],model.params);


xnext = funeval(funfitxy(interp.fspace,interp.Phi,x),interp.fspace,snext);
h     = feval(model.func,'h',ss,xx,[],model.e(repmat(1:K,1,n),:),snext,xnext,model.params);
z     = reshape(model.w'*reshape(h,K,n*p),n,p);

interp.cx = x;
interp.cz = z;
interp.ch = h;

%%%%%%%%%%%%%%%%%%%%%%%%%

options.eqsolver    = 'lmmcp';
%options.eqsolver    = 'path';
options.reesolver   = 'SA';
options.reesolver   = 'mixed';
options.method      = 'resapprox-complete';


tic;
[c,x,z] = recsSolveREE(interp,model,s,x,options);
toc;


% I don't really care about the results now

return


%% plot results


%%%%%%%%%%%%

cont0 = @(sl) funeval(xinit,interp.fspace,[1, sl]);

xt = x;
cont = @(sl) funeval(xt,interp.fspace,[1, sl]);



xxx = 0.5:0.01:5.0


yyy0 = [xxx*0; xxx*0];
yyy = [xxx*0; xxx*0];
for i = 1:length(xxx)
    yyy0(:,i) = cont0(xxx(i));  % perturbation solution
    yyy(:,i) = cont(xxx(i));
end

figure(1)

subplot(211)
plot(xxx,xxx,'black')
hold on
plot(xxx,yyy0(1,:),'green')
plot(xxx,yyy(1,:)) 
title('Consumption')

subplot(212)
plot(xxx,yyy0(2,:),'green') 
hold on
plot(xxx,yyy(2,:))  
title('Welfare')
