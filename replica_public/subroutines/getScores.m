function [Fout,Frawout,Lout] = getScores(in,k)
% Estimate principal components via singular value decomposition (SVD) of
% the data matrix "in" after having standardised the data. The function
% returns the first k scores and the loadings.


% Get dimensions
[T,N] = size(in);

% Standardise data
mu  = repmat(mean(in),T,1);
sdt = repmat(std(in),T,1);
X   = (in - mu)./sdt;

% Perform singular value decomposition
[U,S,V] = svd(X,0);

Fraw = U*S;                  %raw scores (having variances = eigenvalues)
F = U*sqrt(T-1);             %standardized scores
L = V*S/sqrt(T-1);           %loadings

sgn = diag(sign(L(1,:)));    %normalize sign of all factors to be positively
F = F*sgn;                   %correlated with the first observable data series
L = L*sgn;

% Select first k scores
Frawout = Fraw(:,1:k);
Fout = F(:,1:k);
Lout = L(:,1:k);

% % Alternative using eigen-decomposition
% 
% % Covariance
% C = X'*X/(T-1);
% 
% % Eigen-decomposition
% [evec,eval] = eig(C);
% 
% % Sort the eigenvalues in decreasing order
% [eval,idx] = sort(diag(eval),'descend');
% evc = zeros(N,N);
% for i = 1:N
%    evc(:,i) = evec(:,idx(i));
% end
% 
% % Get raw scores
% Fraw = X*evc;
% 
% % Check that eigenvalues map into singular values
% eval_check = diag(S).^2/(T-1);





