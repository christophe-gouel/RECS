%% STO3 Anticipated switch to a public storage policy
% We consider here the transitional dynamics that occurs when a public storage
% policy defending a floor price is announced a few periods before its
% beginning. This is a finite-horizon problem in which the terminal condition is
% not a unique situation, but a state-contingent policy function. This is a
% transition from problem <sto1.html STO1> to <sto2.html STO2>.

%% Enter parameters common to both models
k     = 0.02;
delta = 0.02;
r     = 0.05;
mu    = 10;
alpha = -0.4;
%%
% Shock distribution
Mu      = 1;
sigma   = 0.05;
[e,w]   = qnwnorm(5,Mu,sigma^2);

%% Define approximation space
[interp,s] = recsinterpinit(50,min(e)*0.95,2);

%% Terminal condition: solve for the infinite-horizon behavior under policy
% See <sto2.html STO2> for details.
%
% *Pack model structure*
model2.e       = e;
model2.w       = w;
model2.func    = @sto2model;                               
PF             = 1.02;
Sgbar          = 0.4;
model2.params  = {k,delta,r,mu,alpha,PF,Sgbar};

%%
% *Find a first guess through the perfect foresight solution*
[interp2,x2] = recsFirstGuess(interp,model2,s,1,[0 1 1 0],5);

%% 
% *Solve for rational expectations*
[~,x2] = recsSolveREE(interp2,model2,s,x2);

%% Transition: model without policy in finite horizon
% *Pack model structure*
model1.e       = e;
model1.w       = w;
model1.func    = @sto1model;                               
model1.params  = {alpha,k,delta,r,mu};
%%
% The policy is announced 6 periods before it begins (so $T=7$). In the seventh
% period, the equilibrium is the one defined by the problem with public storage:
T          = 7;
xT         = x2(:,1:3);
options.eqsolver = 'ncpsolve'; % lmmcp does not work in this case
[~,X]  = recsSolveREEFiniteHorizon(interp,model1,s,x2(:,1:3),xT,T,options);

%% Initial period: model without policy under infinite horizon
% In period 0, the policy is not expected by the agents and the prevailing
% equilibrium is the one corresponding to the infinite-horizon competitive
% storage problem (<sto1.html STO1>):
[~,x1] = recsSolveREE(interp,model1,s,X(:,1:3,1));

%% Plot price function
figure
plot(s,[x1(:,3) squeeze(X(:,3,:))])
legend([repmat('Period ',T+1,1),num2str((0:T)')])
legend('Location','NorthEast')
legend('boxoff')
xlabel('Availability')
ylabel('Price')
title('Price function in each period')

%%
%
% <html>
% <hr/>
% </html>
%
% Copyright (C) 2011-2012 Christophe Gouel
%
% Licensed under the Expat license, see <LICENSE.txt>


