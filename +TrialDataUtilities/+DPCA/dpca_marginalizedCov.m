function covs = dpca_marginalizedCov(Xfull, varargin)
import(getPackageImportString);

% default input parameters
options = struct('combinedParams', [],       ...   
                 'lambda',         0,        ...
                 'order',          'yes',    ...
                 'timeSplits',     [],       ...
                 'timeParameter',  [],       ...
                 'notToSplit',     []);

% read input parameters
optionNames = fieldnames(options);
if mod(length(varargin),2) == 1
	error('Please provide propertyName/propertyValue pairs')
end
for pair = reshape(varargin,2,[])    % pair is {propName; propValue}
	if any(strcmp(pair{1}, optionNames))
        options.(pair{1}) = pair{2};
    else
        error('%s is not a recognized parameter name', pair{1})
	end
end

% centering
X = Xfull(:,:);
X = bsxfun(@minus, X, nanmean(X,2));
XfullCen = reshape(X, size(Xfull));

% % total variance
% totalVar = nansum(X(:).^2);

% marginalize
[Xmargs, ~] = dpca_marginalize(XfullCen, 'combinedParams', options.combinedParams, ...
                    'timeSplits', options.timeSplits, ...
                    'timeParameter', options.timeParameter, ...
                    'notToSplit', options.notToSplit, ...
                    'ifFlat', 'yes');

% loop over marginalizations
N = size(Xfull, 1);
covs = cell(length(Xmargs), 1);
for i=1:length(Xmargs)
   mY = reshape(Xmargs{i}, N, []);
   mY = mY(:, all(~isnan(mY), 1));

   covs{i} =  mY*mY.'/size(mY,2);
end
%     
%     if length(options.lambda) == 1
%         thisLambda = options.lambda;
%     else
%         thisLambda = options.lambda(margNums(i));
%     end
%     
%     % regularization
%     if thisLambda ~= 0
%         Xr = [X totalVar * thisLambda * eye(size(X,1))];
%         Xf = [Xmargs{i} zeros(size(X,1))];
%     else
%         Xr = X;
%         Xf = Xmargs{i};
%     end
%     
%     % @djoshea nan fix
%     skip = any(isnan(Xr), 1);
%     XrN = Xr(:, ~skip);
%     XfN = Xf(:, ~skip);
%     
%     % matlab's recommended way
%     % C = Xf/Xr;
%     % [U,~,~] = svd(C*Xr);
%     % U = U(:,1:nc);
%     
%     % fast dirty way
%     
%     % catching possible warning
%     s1 = warning('error','MATLAB:singularMatrix');
%     s2 = warning('error','MATLAB:nearlySingularMatrix');
%     try
%         % fixing C
%         C = XfN*XrN'/(XrN*XrN');
%     catch exception
%         display('Matrix close to singular, using tiny regularization, lambda = 1e-10')
%         thisLambda = 1e-10;
%         Xr = [X totalVar * thisLambda * eye(size(X,1))];
%         Xf = [Xmargs{i} zeros(size(X,1))];
%         
%         % @djoshea nan fix
%         skip = any(isnan(Xr), 1);
%         XrN = Xr(:, ~skip);
%         XfN = Xf(:, ~skip);
%         
%         C = XfN*XrN'/(XrN*XrN');
%     end
%     warning(s1)
%     warning(s2)
%     
%     covs{i} = C;
% end

end