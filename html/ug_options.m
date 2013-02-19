%% Options
%
%% Options structure
% Most RECS functions can have their behavior changed by adjusting
% options. Options take the form of a structure whose possible fields are
% described in the following table. You can set values of these fields using
% either of the two formulations:
%
%  options = struct('field1', values1, 'field2', values2, ...);
%
% or
%
%  options.field1 = values1;
%  options.field2 = values2;
%
% <html>
%   <table border="1">
%     <caption>List of RECS options</caption>
%     <tr>
%       <th> Field name </th> <th> Values </th> <th> Used by Functions </th>
%     </tr>
%     <tr>
%       <td> <tt>accuracy</tt> </td>
%       <td> <tt>1</tt> to check accuracy on the asymptotic distribution (default: <tt>0</tt>). </td>
%       <td>  <a href="matlab:doc('recsSimul')"> <tt>recsSimul</tt> </a> </td>
%     </tr>
%     <tr>
%       <td> <tt>display</tt> </td>
%       <td> <tt>1</tt> (default) to display output (solver iterations or results), <tt>0</tt> to prevent display. </td>
%       <td> <a href="matlab:doc('recsAccuracy')"><tt>recsAccuracy</tt></a>, <a href="matlab:doc('recsmodelinit')"><tt>recsmodelinit</tt></a>, <a href="matlab:doc('recsSimul')"><tt>recsSimul</tt></a>, <a href="matlab:doc('recsSolveREE')"><tt>recsSolveREE</tt></a> </td>
%     </tr>
%     <tr>
%       <td> <tt>eqsolver</tt> </td>
%       <td> Solver for equilibrium equations: '<tt>fsolve</tt>', '<tt>lmmcp</tt>' (default), '<tt>ncpsolve</tt>' or '<tt>path</tt>', see <a href="ug_solvers_eq.html"> Solvers for equilibrium equations</a>. </td>
%       <td>  <a href="matlab:doc('recsFirstGuess')"><tt>recsFirstGuess</tt></a>, <a href="matlab:doc('recsmodelinit')"><tt>recsmodelinit</tt></a>, <a href="matlab:doc('recsSimul')"><tt>recsSimul</tt></a>, <a href="matlab:doc('recsSolveREE')"><tt>recsSolveREE</tt></a>, <a href="matlab:doc('recsSS')"><tt>recsSS</tt></a> </td>
%     </tr>
%     <tr>
%       <td> <tt>eqsolveroptions</tt> </td>
%       <td> Options structure to be passed to eqsolver. </td>
%       <td> <a href="matlab:doc('recsFirstGuess')"><tt>recsFirstGuess</tt></a>, <a href="matlab:doc('recsmodelinit')"><tt>recsmodelinit</tt></a>, <a href="matlab:doc('recsSimul')"><tt>recsSimul</tt></a>, <a href="matlab:doc('recsSolveREE')"><tt>recsSolveREE</tt></a>, <a href="matlab:doc('recsSS')"><tt>recsSS</tt></a> </td>
%     </tr>
%     <tr>
%       <td> <tt>extrapolate</tt> </td>
%       <td> <tt>1</tt> (default) if extrapolation is allowed outside the interpolation space or <tt>0</tt> to forbid it. </td>
%       <td> <a href="matlab:doc('recsSimul')"><tt>recsSimul</tt></a>, <a href="matlab:doc('recsSolveREE')"><tt>recsSolveREE</tt></a> </td>
%     </tr>
%     <tr>
%       <td> <tt>funapprox</tt> </td>
%       <td> Functional approximated to find the rational expectations solution:
%         <ul>
%           <li> '<tt>expapprox</tt>': expectations (variables <tt>z</tt>), similar to the parameterized expectations approach of den Haan and Marcet (1990).</li>
%           <li> '<tt>expfunapprox</tt>': expectations function (function <tt>h</tt>).</li>
%           <li> '<tt>resapprox-simple</tt>' or '<tt>resapprox-complete</tt>' (default): response variables (variables <tt>x</tt>).</li>
%         </ul>
%       </td>
%       <td> <a href="matlab:doc('recsSimul')"><tt>recsSimul</tt></a>, <a href="matlab:doc('recsSolveREE')"><tt>recsSolveREE</tt></a> </td>
%     </tr>
%     <tr>
%       <td> <tt>functional</tt> </td>
%       <td> <tt>1</tt> if the equilibrium equations are a functional equation problem (default: <tt>0</tt>), used to solve optimal discretionary policy. </td>
%       <td> <a href="matlab:doc('recsSimul')"><tt>recsSimul</tt></a>, <a href="matlab:doc('recsSolveREE')"><tt>recsSolveREE</tt></a> </td>
%     </tr>
%     <tr>
%       <td> <tt>loop_over_s</tt> </td>
%       <td> <tt>0</tt> (default) to solve all grid points at once or <tt>1</tt> to loop over each grid points. </td>
%       <td> <a href="matlab:doc('recsSolveREE')"><tt>recsSolveREE</tt></a> </td>
%     </tr>
%     <tr>
%       <td> <tt>Python</tt> </td>
%       <td> 1 to call Python directly instead of the executable file (default: 0, only for Windows and for developement).</td>
%       <td> <a href="matlab:doc('recsmodelinit')"><tt>recsmodelinit</tt></a> </td>
%     </tr>
%     <tr>
%       <td> <tt>reemethod</tt> </td>
%       <td> Define the method to solve for rational expextations:
%            <ul>
%              <li> '<tt>iter</tt>' (default): time iteration.</li>
%              <li> '<tt>iter-newton</tt>': time iteration with Newton step.</li>
%              <li> '<tt>1-step</tt>': rational expectations and equilibrium equations are solved in one step (as suggested in Judd, 1992).</li>
%            </ul>
%       </td>
%       <td> <a href="matlab:doc('recsSolveREE')"><tt>recsSolveREE</tt></a> </td>
%     </tr>
%     <tr>
%       <td> <tt>reesolver</tt> </td>
%       <td> Solver for the rational expectations equilibrium, expected values depend on <tt>reemethod</tt>:
%            <ul>
%              <li> if <tt>iter</tt>: '<tt>kinsol</tt>', '<tt>krylov</tt>', '<tt>mixed</tt>', or '<tt>SA</tt>' (default).</li>
%              <li> if <tt>iter-newton</tt>: '<tt>fsolve</tt>', '<tt>lmmcp</tt>', '<tt>ncpsolve</tt>' or '<tt>path</tt>', see <a href="ug_solvers_eq.html">Solvers for equilibrium equations</a>.</li>
%              <li> if <tt>1-step</tt>: irrelevant, the solver is determined by <tt>eqsolver</tt>.</li>
%            </ul>
%       </td>
%       <td> <a href="matlab:doc('recsSolveREE')"><tt>recsSolveREE</tt></a> </td>
%     </tr>
%     <tr>
%       <td> <tt>reesolveroptions</tt> </td>
%       <td> Options structure to be passed to reesolver. </td>
%       <td> <a href="matlab:doc('recsSolveREE')"><tt>recsSolveREE</tt></a> </td>
%     </tr>
%     <tr>
%       <td> <tt>simulmethod</tt> </td>
%       <td> Method of simulation: '<tt>interpolation</tt>' (default) or '<tt>solve</tt>', see <a href="simulate.html">Simulate the model</a>. </td>
%       <td>  <a href="matlab:doc('recsSimul')"><tt>recsSimul</tt></a> </td>
%     </tr>
%     <tr>
%       <td> <tt>stat</tt> </td>
%       <td> <tt>1</tt> to ouput summary statistics from the simulation (default: <tt>0</tt>), see <a href="simulate.html">Simulate the model</a>. </td>
%       <td> <a href="matlab:doc('recsSimul')"><tt>recsSimul</tt></a> </td>
%     </tr>
%     <tr>
%       <td> <tt>useapprox</tt> </td>
%       <td> (default: <tt>1</tt>) behaviour dependent of the chosen function to
%             approximate. If <tt>0</tt> and funapprox is '<tt>expapprox</tt>' then next-period
%             responses are calculated by equations solve and not just
%             interpolated. If <tt>1</tt> and funapprox is '<tt>resapprox</tt>', the guess of
%             response variables is found with the new approximation structure.
%       </td>
%       <td> <a href="matlab:doc('recsSimul')"><tt>recsSimul</tt></a>, <a href="matlab:doc('recsSolveREE')"><tt>recsSolveREE</tt></a> </td>
%     </tr>
%   </table>
% </html>
%
%% References
% <http://www.jstor.org/stable/1391746 den Haan, W. J. and Marcet, A.
% (1990). Solving the Stochastic Growth Model by Parameterizing
% Expectations. _Journal of Business & Economic Statistics_, 8(1), 31-34.>
%
% <http://dx.doi.org/10.1016/0022-0531(92)90061-L Judd, K. (1992). Projection
% Methods for Solving Aggregate Growth Models. _Journal of Economic Theory_,
% 58(2), 410-452.>
