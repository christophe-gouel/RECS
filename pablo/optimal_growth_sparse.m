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


% DEFINE DEFINITION SET

smin  = [    (1-sigma1*1.5)  0.5 ];     % minimum state
smax  = [    (1+sigma1*1.5)  2 ];     % maximum state

%% Define a sparse grid
% d: number of dimensions
% l: maximum level of subgrids
% gridtype : type of sparse grid (any type accepted by the sparse toolbox)
l = 4;
d = 2;
gridtype = 'Clenshaw-Curtis';
%gridtype = 'Maximum';
%gridtype = 'Gauss-Patterson';
%gridtype = 'NoBoundary';
%gridtype = 'Chebyshev';
interp = create_sp_grid(d, l, gridtype, smin,smax);
s = interp.s;
n = size(s,1);
% suppress unavoidable warnings
warning('off','MATLAB:spinterp:outOfRange');
warning('off','MATLAB:spinterp:insufficientDepth');


%% Define a functional space with compecon toolbox
%interp = struct;
%order     = [ 5 50  ];                             % degree of approximation
%interp.fspace = fundefn('cheb',order,smin,smax);       % function space
%interp.fspace = fundef( {'cheb',5,smin(1),smax(1)}, {'spli',[smin(2),smax(2)],50});       % function space
%snodes = funnode(interp.fspace);                   % state collocation nodes
%s=gridmake(snodes);
%n             = size(s,1);
%interp.Phi    = funbasx(interp.fspace);
%interp.grid = s;
%interp.smin = smin;
%interp.smax = smax;
%interp.type = 'compecon';


%% initial guess

dsol = feval(model.func,'solution'); % deterministic solution
[LB,UB] = feval(model.func,'b',s); % bounds

s_ss = repmat(dsol.s_ss,n,1);
x_ss = repmat(dsol.x_ss,n,1);
xinit = x_ss + (s - s_ss) * dsol.Z';

xinit = max(min(xinit,UB-0.001),LB+0.001);
%xinit = max(min(xinit,UB),LB);

%% options

options.eqsolver    = 'lmmcp';
%options.eqsolver    = 'path';
options.reesolver   = 'SA';

options.reesolveroptions.showiters   = 1;

options.method      = 'resapprox-complete';




tic;
[c,x,z] = pabSolveREE(interp,model,s,xinit,options);
toc;

return 
tic;
interp.cx = fit_coefficients(interp,xinit);
%interp.cx = xinit;
[c,x,z] = recsSolveREE(interp,model,s,xinit,options);
toc;


% I don't really care about the results now

return


%% plot results


%%%%%%%%%%%%


c0 = fit_coefficients(interp, xinit);

cont0 = @(sl) eval_function(interp,c0,[1, sl]);

xt = x;


cont = @(sl) eval_function(interp,c,[1, sl]);


xxx = 0.5:0.01:2.0;


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
