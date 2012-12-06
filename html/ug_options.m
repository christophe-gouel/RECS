%% Options
%
%
% <html>
%   <table border="1">
%     <caption>List of RECS options</caption>
%     <tr>
%       <th> <tt>Field name </th> <th> Values </th> <th> Used by Functions </th>
%     </tr>
%     <tr>
%       <td> <tt>accuracy</tt> </td>
%       <td> 1 to check accuracy on the asymptotic distribution (default: 0) </td>
%       <td>  <a href="matlab:doc('recsSimul')"> <tt>recsSimul</tt> </a> </td>
%     </tr>
%     <tr>
%       <td> <tt>display</tt> </td>
%       <td> 1 to display output (solver iterations or results) (default: 1) </td>
%       <td> <a href="matlab:doc('recsmodelinit')"> <tt>recsmodelinit</tt> </a>, <a href="matlab:doc('recsSolveREE')"> <tt>recsSolveREE</tt> </a> </td>
%     </tr>
%     <tr>
%       <td> <tt>eqsolver</tt> </td>
%       <td> '<tt>fsolve</tt>', '<tt>lmmcp</tt>' (default), '<tt>ncpsolve</tt>' or '<tt>path</tt>' </td>
%       <td>  <a href="matlab:doc('recsFirstGuess')"> <tt>recsFirstGuess</tt> </a>, <a href="matlab:doc('recsmodelinit')"> <tt>recsmodelinit</tt> </a>,<a href="matlab:doc('recsSimul')"> <tt>recsSimul</tt> </a>, <a href="matlab:doc('recsSolveREE')"> <tt>recsSolveREE</tt> </a>, <a href="matlab:doc('recsSS')"> <tt>recsSS</tt> </a> </td>
%     </tr>
%     <tr>
%       <td> <tt>eqsolveroptions</tt> </td>
%       <td> options structure to be passed to eqsolver </td>
%       <td> <a href="matlab:doc('recsFirstGuess')"> <tt>recsFirstGuess</tt> </a>, <a href="matlab:doc('recsmodelinit')"> <tt>recsmodelinit</tt> </a>,<a href="matlab:doc('recsSimul')"> <tt>recsSimul</tt> </a>, <a href="matlab:doc('recsSolveREE')"> <tt>recsSolveREE</tt> </a>, <a href="matlab:doc('recsSS')"> <tt>recsSS</tt> </a> </td>
%     </tr>
%     <tr>
%       <td> <tt>extrapolate</tt> </td>
%       <td> 1 if extrapolation is allowed outside the interpolation space or 0 to forbid it (default: 1) </td> 
%       <td> <a href="matlab:doc('recsSimul')"> <tt>recsSimul</tt> </a>, <a href="matlab:doc('recsSolveREE')"> <tt>recsSolveREE</tt> </a> </td>
%     </tr>
%     <tr>
%       <td> <tt>funapprox</tt> </td> 
%       <td> '<tt>expapprox</tt>', '<tt>expfunapprox</tt>', '<tt>resapprox-simple</tt>' or '<tt>resapprox-complete</tt>' (default) </td>
%       <td> <a href="matlab:doc('recsSimul')"> <tt>recsSimul</tt> </a>, <a href="matlab:doc('recsSolveREE')"> <tt>recsSolveREE</tt> </a> </td>
%     </tr>
%     <tr>
%       <td> <tt>functional</tt> </td> 
%       <td> 1 if the equilibrium equations are a functional equation problem (default: 0) </td>
%       <td> <a href="matlab:doc('recsSimul')"> <tt>recsSimul</tt> </a>, <a href="matlab:doc('recsSolveREE')"> <tt>recsSolveREE</tt> </a> </td>
%     </tr>
%     <tr>
%       <td> <tt>loop_over_s</tt> </td>
%       <td> 0 (default) to solve all grid points at once or 1 to loop over each grid points </td>
%       <td> <a href="matlab:doc('recsSolveREE')"> <tt>recsSolveREE</tt> </a> </td>
%     </tr>
%     <tr>
%       <td> <tt>reemethod</tt> </td> 
%       <td> '<tt>iter</tt>' (default), '<tt>iter-newton</tt>' or '<tt>1-step</tt>' </td>
%       <td> <a href="matlab:doc('recsSolveREE')"> <tt>recsSolveREE</tt> </a> </td>
%     </tr>
%     <tr>
%       <td> <tt>reesolver</tt> </td>
%       <td> '<tt>krylov</tt>', '<tt>mixed</tt>', '<tt>SA</tt>' (default) or '<tt>fsolve</tt>' </td>
%       <td> <a href="matlab:doc('recsSolveREE')"> <tt>recsSolveREE</tt> </a>  </td>
%     </tr>
%     <tr>
%       <td> <tt>reesolveroptions</tt> </td> 
%       <td> options structure to be passed to reesolver </td> 
%       <td> <a href="matlab:doc('recsSolveREE')"> <tt>recsSolveREE</tt> </a> </td>
%     </tr>
%     <tr>
%       <td> <tt>simulmethod</tt> </td>
%       <td> '<tt>interpolation</tt>' (default) or '<tt>solve</tt>' </td>
%       <td>  <a href="matlab:doc('recsSimul')"> <tt>recsSimul</tt> </a> </td>
%     </tr>
%     <tr>
%       <td> <tt>stat</tt> </td>
%       <td> 1 to ouput summary statistics from the simulation (default: 0) </td>
%       <td> <a href="matlab:doc('recsSimul')"> <tt>recsSimul</tt> </a> </td>
%     </tr>
%     <tr>
%       <td> <tt>useapprox</tt> </td>
%       <td> (default: 1) behaviour dependent of the chosen function to
%             approximate. If 0 and funapprox is '<tt>expapprox</tt>' then next-period
%             responses are calculated by equations solve and not just
%             interpolated. If 1 and funapprox is '<tt>resapprox</tt>', the guess of
%             response variables is found with the new approximation structure
%       </td>
%       <td> <a href="matlab:doc('recsSimul')"> <tt>recsSimul</tt> </a>, <a href="matlab:doc('recsSolveREE')"> <tt>recsSolveREE</tt> </a> </td>
%     </tr>
%   </table>
% </html>