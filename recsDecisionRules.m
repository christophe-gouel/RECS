function recsDecisionRules(model,interp,states,space,s0)

%% Initialization
if nargin<3 || isempty(states), states = 1:interp.fspace.d; end
if nargin<4 || isempty(space)
  range = [interp.fspace.a(states)' interp.fspace.b(states)'];
  n     = 100*ones(1,length(states));
else
  range = space(1:2,:);
  n     = space(3,:);
end
if nargin<5 || isempty(s0), s0 = model.sss; end

d = size(s0,2);

for i=1:length(states)
  s              = repmat(s0,1,n(i));
  s              = reshape(s',d,[])';
  s(:,states(i)) = repmat(linspace(range(i,1),range(i,2),n(i)),...
                          1,...
                          size(s0,1));

  [~,x] = recsSimul(model,interp,s,0,[]);

  m = size(x,2);
  x = reshape(x,n(i),[],m);
  x = permute(x,[1 3 2]);
  figure
  for j=1:m
    subplot(ceil((m)/ceil(sqrt(m))),ceil(sqrt(m)),j)
    plot(s(1:n(i),states(i)),squeeze(x(:,j,:)))
  end
end
