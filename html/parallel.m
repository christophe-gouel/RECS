%% Parallel computing
% If you have access to MATLAB Parallel Computing Toolbox, you can use parallel
% computing to speed up RECS. Parallel computing may matter a lot because models
% are solved on grids and, with most solution methods, the solution for each
% grid point can be solved independently from other grid points. It mostly
% matters in three steps: (i) the calculation of first guess through perfect
% foresight, (ii) the solution of stochastic rational expectations problems, and
% (iii) simulating models using equations solvers (rather than interpolation).
%
%% Enabling parallel computing
% For MATLAB versions before R2013b, parallel computation is enabled for MATLAB
% Parallel Computing Toolbox by the following call
%
%  matlabpool
%
% After MATLAB R2013b, by default, parallel computation automatically starts
% when required or can be started using
%
%  parpool
%
%
%% Parallel computing for computing first guess through perfect foresight
% Solving for a first guess through perfect foresight can take as much time as
% solving for the stochastic rational expectations problem because a
% perfect-foresight solution has to be solved for each grid point
% independently.
%
% As soon as MATLAB parallel computing features are enabled, parallel computing
% will be used automatically by |recsFirstGuess| by using each MATLAB worker to
% solve for one grid point.
%
%% Parallel computing for solving for stochastic rational expectations equilibrium
% Parallel computing can be used with |recsSolveREE| and
% |recsSolveREEFiniteHorizon|. By default, these functions solve for the
% rational expectations equilibrium by solving for all grid points at once the
% equilibrium equations for given expectations functions and by iterating
% through time. Grid points are solved all at once because it tends to be faster
% to solve a large system of equations than to solve many small systems. To
% benefit from parallel computing, the large system of equations has to be
% broken down in smaller systems that can be fed to the MATLAB workers. This is
% governed by the option |loop_over_s|.
%
% By default |loop_over_s| is equal to 0, which implies that all grid points are
% solved at once; parallel computing will not be used in this case. Setting
% |loop_over_s| equal to 1 implies that each grid point is solved
% separately. This option can benefit from parallel computing, but it usually
% breaks down the problem in too many small systems of equations. For
% |loop_over_s|, any integer value |n| different from 0 and 1 implies that
% equations will be solved by |n| blocks of grid points. To keep each MATLAB
% worker occupied, a good rule is to set |n| equal to a multiple of the number
% of MATLAB workers.
%
% It should also be noted that solving smaller systems of equations is easier
% than solving at once for all grid points (but if it may be slower). So
% increasing |loop_over_s| is, independently from parallel computing, useful to
% facilitate the solution process is the problem is difficult to solve.
%
% In the case of |recsSolveREE|, parallel computing is not compatible with the
% option field |reemethod| equal to |1-step|.
%
%% Simulating models using equations solvers
% There are two techniques to simulate models once the rational expectations
% solution has been found (see <simulate.html Simulate the model>). The most
% precise and slowest technique simulates the model by solving equilibrium
% equations using approximated expectations (option field |simulmethod| set to
% |solve|). If several trajectories are simulated at the same time, the approach
% is similar to the one described above for |recsSolveREE|: by default,
% equilibrium equations for all trajectories are solved at once. To benefit from
% parallel computing, one needs to do as above: set |loop_over_s| to a non-zero
% integer value.
