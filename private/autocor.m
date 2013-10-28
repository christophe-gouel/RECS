function acor = autocor(x)
% AUTOCOR Computes the sample autocorrelation of a matrix of time-series x

maxlag = 5;
[N,d]  = size(x);
acov   = zeros(maxlag+1,d);
for n=0:maxlag
  acov(n+1,:) = mean(x(1+n:N,:).*x(1:N-n,:),1)-mean(x(1+n:N,:),1).*mean(x(1:N-n,:),1);
end

acor                  = (acov(2:maxlag+1,:)./acov(ones(maxlag,1),:))';

% If variance is smaller than precision, correlation is equal to 1
testprecision         = eps(mean(x))>acov(1,:);
acor(testprecision,:) = 1;
