classdef KalmanFilterDecoder < handle
    % Class for the Kalman Filter Decoder
    % translated from https://github.com/KordingLab/Neural_Decoding
    % for manuscript https://arxiv.org/abs/1708.00909
    % Machine learning for neural decoding
    % Joshua I. Glaser, Ari S. Benjamin, Raeed H. Chowdhury, Matthew G. Perich, Lee E. Miller, Konrad P. Kording
    % specifically, https://github.com/KordingLab/Neural_Decoding/blob/ed8e8f92027f92c18ebfd877ef28333da2e9646c/Neural_Decoding/decoders.py#L201
    
    % Parameters
    % -----------
    % C - float, optional, default 1
    % This parameter scales the noise matrix associated with the transition in kinematic states.
    % It effectively allows changing the weight of the new neural evidence in the current update.
    
    % Our implementation of the Kalman filter for neural decoding is based on that of Wu et al 2003 (https://papers.nips.cc/paper/2178-neural-decoding-of-cursor-motion-using-a-kalman-filter.pdf)
    %   with the exception of the addition of the parameter C.
    %   The original implementation has previously been coded in Matlab by Dan Morris (http://dmorris.net/projects/neural_decoding.html#code)
    
    properties
      C (1, 1) double
      
      A (:, :) double % transition matrix
      W (:, :) double % covariance of transition matrix
      H (:, :) double % measurment matrix
      Q (:, :) double % covariance of measurement matrix
    end
    
    methods
      function kf = KalmanFilterDecoder(args)
        arguments
          args.C (1, 1) = 1;
        end
        
        kf.C = args.C;
      end
      
      function fit(kf, X_kf_train, y_train)
        % Train Kalman Filter Decoder
        
        % Parameters
        % ----------
        % X_kf_train: numpy 2d array of shape [n_samples(i.e. timebins) , n_neurons]
        % This is the neural data in Kalman filter format.
        % See example file for an example of how to format the neural data correctly
        
        % y_train: numpy 2d array of shape [n_samples(i.e. timebins), n_outputs]
        % This is the outputs that are being predicted
        
        % First we'll rename and reformat the variables to be in a more standard kalman filter nomenclature (specifically that from Wu et al, 2003):
        % xs are the state (here, the variable we're predicting, i.e. y_train)
        % zs are the observed variable (neural data here, i.e. X_kf_train)
        arguments
          kf
          X_kf_train (:, :)
          y_train (:, :)
        end
        
        X = y_train'; % neurons x time
        Z = X_kf_train'; % outputs x time
        
        % number of time bins
        nt = size(X, 2);
        
        % Calculate the transition matrix (from x_t to x_t+1) using least-squares, and compute its covariance
        % In our case, this is the transition from one kinematic state to the next
        X2 = X(:,2:end);
        X1 = X(:,1:end-1);
        A = (X2*X1') / (X1*X1');  %#ok<*PROPLC> % Transition matrix
        W = (X2-A*X1)*(X2-A*X1)'/(nt-1) / kf.C; % Covariance of transition matrix.
        % Note we divide by nt-1 since only nt-1 points were used in the computation (that's the length of X1 and X2). We also introduce the extra parameter C here.
        
        % Calculate the measurement matrix (from x_t to z_t) using least-squares, and compute its covariance
        % In our case, this is the transformation from kinematics to spikes
        H = (Z*X') / (X*X'); % Measurement matrix
        Q = ((Z - H*X)*((Z - H*X)')) / nt; % Covariance of measurement matrix
        kf.A = A;
        kf.W = W;
        kf.H = H;
        kf.Q = Q;
      end
      
      function y_test_predicted = predict(kf, X_kf_test, initial_state)
        % Predict outcomes using trained Kalman Filter Decoder
        
        % Parameters
        % ----------
        % X_kf_test: numpy 2d array of shape [n_samples(i.e. timebins) , n_neurons]
        % This is the neural data in Kalman filter format.
        
        % y_test: numpy 2d array of shape [n_samples(i.e. timebins),n_outputs]
        % The actual outputs
        % This parameter is necesary for the Kalman filter (unlike other decoders)
        % because the first value is nececessary for initialization
        
        % Returns
        % -------
        % y_test_predicted: numpy 2d array of shape [n_samples(i.e. timebins),n_outputs]
        % The predicted outputs
        arguments
          kf
          X_kf_test (:, :)
          initial_state (1, :)
        end
        
        % #Extract parameters
        A = kf.A;
        W = kf.W;
        H = kf.H;
        Q = kf.Q;
        
        % First we'll rename and reformat the variables to be in a more standard kalman filter nomenclature (specifically that from Wu et al):
        % xs are the state (here, the variable we're predicting, i.e. y_train)
        % zs are the observed variable (neural data here, i.e. X_kf_train)
        
        %X = y_test';
        Z = X_kf_test';
        num_states = size(A, 1); % Dimensionality of the state
        num_time = size(Z, 2);
        
        % Initializations
        states = nan(num_states, num_time); % Keep track of states over time (states is what will be returned as y_test_predicted)
        state = initial_state; % Initial state
        states(:, 1) = initial_state;
        
        % Get predicted state for every time bin
        for t = 1:num_time
          % Do first part of state update - based on transition matrix
          P_m = A*P*A' + W; % num_states x num_states
          state_m = A*state;
          
          % Do second part of state update - based on measurement matrix
          K = (P_m*H') / (H*P_m*H'+Q); % Calculate Kalman gain
          P = (eye(num_states)-K*H)*P_m; % num_states x num_states
          state = state_m + K*(Z(:,t+1)-H*state_m);
          states(:,t+1) = state; % Record state at the timestep
        end
        
        y_test_predicted = states';
      end
    end
end



