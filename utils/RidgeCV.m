classdef RidgeCV < handle
    properties
        lambda_search (:, 1) = logspace(-8,-1,20);
        shared_lambda (1, 1) logical = true;
        kfold (1, 1) {mustBeInteger} = 5;
        standardize_cols_individually (1, 1) logical = false;
    end

    properties
        lambda (1, :)

        beta (:, :) % nFeatures (cols in X) x nObservations (cols in y)
        bias (1, :) % 1 x nObservations (cols in y)
    end

    methods
        function fit(rcv, X, Y)
            arguments
                rcv
                X (:, :)
                Y (:, :)
            end

            assert(size(X,1) == size(Y, 1), 'Number of observations (rows) unequal');

            % standardize all columns collectively 
            if rcv.standardize_cols_individually
                transform_scale = std(X, 0, 1, "omitnan");
            else
                transform_scale = std(X, 0, "all", "omitnan");
            end
        
            X = X ./ transform_scale;

            ID = size(X, 2);
            OD = size(Y, 2);
            nLambda = numel(rcv.lambda_search);

            % strip observations with missing valid
            valid_mask = all(~isnan(X), 2) & all(~isnan(Y), 2);
            X = X(valid_mask, :);
            Y = Y(valid_mask, :);

            mse = nan(nLambda, OD);
            
            for od = 1:OD
                y = Y(:, od);
                cv_mdl = fitrlinear(X', y', ObservationsIn='columns', KFold=rcv.kfold, ...
                    Lambda=rcv.lambda_search, Learner='leastsquares', ...
                    Solver='lbfgs',Regularization='ridge');

                mse(:, od) = kfoldLoss(cv_mdl);
            end

            % pick lambda minimizing cv MSE
            if rcv.shared_lambda
                mse = sum(mse, 2);
                [~, ind_lambda] = min(mse);
                rcv.lambda = repmat(rcv.lambda_search(ind_lambda), 1, OD);
            else
                rcv.lambda = nan(1, OD);
                for od = 1:OD
                    [~, ind_lambda] = min(mse(:, od));
                    rcv.lambda(od) = rcv.lambda_search(ind_lambda);
                end
            end

            % fit final model
            beta_unscaled = nan(ID, OD, like=X);
            bias = nan(1, OD, like=X); %#ok<*PROPLC>
            for od = 1:OD
                y = Y(:, od);
                mdl = fitrlinear(X', y', ObservationsIn='columns', ...
                    Lambda=rcv.lambda(od), Learner='leastsquares', ...
                    Solver='lbfgs',Regularization='ridge');

                beta_unscaled(:, od) = mdl.Beta;
                bias(:, od) = mdl.Bias;
            end

            % apply scaling to beta to simplify transform
            rcv.beta = beta_unscaled ./ (transform_scale');
            rcv.bias = bias;
        end

        function Y = predict(rcv, X)
            arguments
                rcv
                X (:, :)
            end

            % X (N x ID) * beta (ID x OD) + bias (1, OD) == Y (N x OD)
            Y = X * rcv.beta + rcv.bias;
        end
    end

end