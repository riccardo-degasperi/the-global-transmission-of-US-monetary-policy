function Yf = SSAgetBands(Y)
% Get confidence regions out of SSA (time-series) draws

    % Get dimensions
    [h,N,draws] = size(Y);

    % Initialise
    Yf = nan(h,N,5);

    % Specify intervals
    c1 = 0.9;         %90% confidence bands
    c2 = 0.68;        %68% confidence bands


    for ii = 1:N
      for jj = 1:h

          Yf(jj,ii,1) = quantile(Y(jj,ii,:),(1-c1)/2);
          Yf(jj,ii,2) = quantile(Y(jj,ii,:),(1-c2)/2);
          Yf(jj,ii,3) = quantile(Y(jj,ii,:),0.5);
          Yf(jj,ii,4) = quantile(Y(jj,ii,:),1-(1-c2)/2);
          Yf(jj,ii,5) = quantile(Y(jj,ii,:),1-(1-c1)/2);

      end
    end

end