function IRFulm = IRFbands(IRF,N,h)


% Initialise container of IRF median, upper and lower bounds
IRFulm = cell(N,N);

c1 = 0.9;   %90% confidence bands
c2 = 0.68;  %68% confidence bands

% Compute its elements
for i=1:N
   for ii=1:N
      for j=1:h+1

      if isempty(IRF{i,ii})
          %skip
      else
          
          % Lower bound 90%
          IRFulm{i,ii}(1,j) = quantile(IRF{i,ii}(j,:),(1-c1)/2);

          % Lower bound 68%
          IRFulm{i,ii}(2,j) = quantile(IRF{i,ii}(j,:),(1-c2)/2);

          % Median
          IRFulm{i,ii}(3,j) = quantile(IRF{i,ii}(j,:),0.5);

          % Upper bound 68%
          IRFulm{i,ii}(4,j) = quantile(IRF{i,ii}(j,:),1-(1-c2)/2);

          % Upper bound 90%
          IRFulm{i,ii}(5,j) = quantile(IRF{i,ii}(j,:),1-(1-c1)/2);
      
      end
      
      end
   end
end


