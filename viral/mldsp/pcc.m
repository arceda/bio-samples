%PEARSON CORRELATION COEFICIENTS DISTANCE IN POCTAVE, THE VECTOR ARE THE MAGNITUD SPECTRUM OF FFT

	



%X = [1460, 517.201, 230.163, 453.649, 266.169, 267.257; 1340, 569.351, 219.907, 473.615, 239.587, 229.557; 1462, 622.617, 324.276, 503.927, 214.432,223.652];
X = [1460, 517.201, 230.163, 453.649, 266.169, 267.257; 1340, 569.351, 219.907, 473.615, 239.587, 229.557; 1462, 622.617, 324.276, 503.927, 214.432,223.652; 1994, 672.012, 456.685, 412.211, 219.068, 131.52];
X

n = size(X,1); %rows
D   = repmat(NaN,n-1,1); % preallocate
D
idx = [0 0];             % initialize

for i = 1:n-1 % repeat for each row except the last
	 blk    = numel(i+1:n);              % get size of this block of distances
	 r      = repmat(X(i,:),blk,1);      % extract this row, replicate
	 R      = X(i+1:n,:);                % extract remaining rows
	 idx    = idx(end)+1:idx(end)+blk;   % get index to this block of distances
 
	disp('r and R');
	r
	R
	sqrt( sum( (r-R).^2 ,2) )

   	D(idx) = sqrt( sum( (r-R).^2 ,2) ); % sum by row prodces a column
 
end

D
idx

x = D;

%% re wrap matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x     = x(:);      % make sure it's a column vector
sizeX = length(x); % size of tridiagonal

% Get dimensions of symmetric matrix:
% There are [n(n-1)/2] values in the lower tridiagonal, so solve for n:
dim = 0.5*(1+sqrt(8*sizeX+1));


xDis     = zeros(dim,dim);              % preallocate results matrix
ai       = logical(tril(ones(dim),-1)); % get indices for lower diagonal
xDis(ai) = x;                           % fill lower tridiagonal
xDis     = xDis + xDis';                % fill upper tridiagonal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xDis


D = corrcoef(X')
D = ( 1 - (corrcoef(X')) )/2