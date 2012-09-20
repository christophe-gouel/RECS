function B = Bounds(func,s0,params,output,ix)
%% BOUNDS Allows differentiation of bounds

Big = 1E20;
[LBx,UBx]     = func('b',s0,[],[],[],[],[],params);

if output==1
  B = LBx(:,ix(:,1));
else
  B = UBx(:,ix(:,2));
end

B(isinf(B)) = sign(B(isinf(B)))*Big;
B           = B(:);
