% Get data from user.data file
A = importdata('C:\Users\RITESH\Desktop\ee239 final\ml-100k\u.data');
user_id = A(:, 1);
movie_id = A(:, 2);
% rating = A(:, 3);
% 
% % ------------------------------------------------------
% % Part 1
% % ------------------------------------------------------
% 
% % Build matrix R and w (weights matrix)
% R = zeros(943, 1682);
% w = zeros(943, 1682);
% for i=1:100000
%     R(user_id(i), movie_id(i)) = rating(i);
%     w(user_id(i), movie_id(i)) = 1;
% end
% 
% % Three different values of k
% k = [10, 50, 100];
% 
% % Factorize using wnmfrule
% [U11,V11,numIter,tElapsed,finalResidual11]=mywnmfrule(R, k(1), w);
% [U12,V12,numIter,tElapsed,finalResidual12]=mywnmfrule(R, k(2), w);
% [U13,V13,numIter,tElapsed,finalResidual13]=mywnmfrule(R, k(3), w);
% 
% least_squared_error1(1,1) = finalResidual11; 
% least_squared_error1(2,1) = finalResidual12; 
% least_squared_error1(3,1) = finalResidual13; 
% 
% least_squared_error1
% 
% % Final Residual Errors
% %   247.3379
% %   173.9582
% %   131.1707
%   
% % -----------------------------------------------
% % Part 2
% % -----------------------------------------------
% % Swap R and w
% temp = R;
% R2 = w;
% w2 = temp;
% 
% % Three different values of k
% k = [10, 50, 100];
% 
% % Factorize using wnmfrule
% [U21,V21,numIter,tElapsed,finalResidual21]=mywnmfrule(R2,k(1), w2);
% [U22,V22,numIter,tElapsed,finalResidual22]=mywnmfrule(R2,k(2), w2);
% [U23,V23,numIter,tElapsed,finalResidual23]=mywnmfrule(R2,k(3), w2);
% 
% least_squared_error2(1,1) = finalResidual21; 
% least_squared_error2(2,1) = finalResidual22; 
% least_squared_error2(3,1) = finalResidual23;
% 
% least_squared_error2
% 
% % Final Residual Errors
% %    18.2706
% %    26.1834
% %    23.0279
% 
% % -----------------------------------------------
% % Part 3
% % -----------------------------------------------
% random_indices = randperm(100000);
% % Training rating and weight matrices
% R_Actual31 = zeros(943, 1682);
% w31 = zeros(943, 1682);
% % Matrix to store predicted values
% R_Predicted31= zeros(size(R_Actual31, 1), size(R_Actual31, 2), 3); % 3 copies for k=10, 50, 100
% % Matrix to store average error
% average_error31 = zeros(10, 3); %10 folds, 3 values of k
% 
% % k can take the values 10, 50, 100
% k=[10, 50, 100];
% 
% % Training and Testing Data as per 10 - fold Cross Validation
% for j = 1:10
%     
%     % Create R_Actual31 and w31 as the training ratings and weights matrix
%     % from 1/10th to (j-1)th part set as rating, and 1
%     for i = 1:(j-1)*10000
%         R_Actual31(user_id(random_indices(i)), movie_id(random_indices(i))) = rating(random_indices(i));
%         w31(user_id(random_indices(i)), movie_id(random_indices(i))) = 1;
%     end
%     
%     %jth block will be full of zeros
%     
%     % from jth part to the end set as rating, and 1
%     for i=(j*10000 + 1):100000
%       R_Actual31(user_id(random_indices(i)), movie_id(random_indices(i))) = rating(random_indices(i));  
%       w31(user_id(random_indices(i)), movie_id(random_indices(i))) = 1;
%     end
%     
%     % Perform factorization, calculate error, and store predicted values
%     for m = 1:3 % For 3 values of k
%         % Perform factorization
%         [U3,V3,numIter,tElapsed,finalResidual]=mywnmfrule(R_Actual31,k(m), w31);
% 
%         %Predicted values
%         product1 = U3*V3;
%         
%         % calculate the sum of absolute error and predicted matrix storage
%         sum_of_errors = 0;
%         for i=((j-1)*10000 + 1):j*10000
%             
%             % Store indices as row3 and col3
%             row3 = user_id(random_indices(i));
%             col3 =  movie_id(random_indices(i));
%             
%             % Compute absolute error and add it to sum
%             sum_of_errors = sum_of_errors + abs(product1(row3, col3) - R(row3, col3));
%             
%             % Store the predicted value
%             R_Predicted31(row3, col3,m) = product1(row3, col3);
%         end
%         
%         % determine the average absolute error 
%         average_error31(j, m) = sum_of_errors / 10000
%     end
% end
% 
% % Now print the average amongst 10 tests for each of the k values
% for i=1:3
%     sprintf('For k = %d value', k(i))
%     
%     %Average error of 10 folds
%     sprintf('Average error = %f', mean(average_error31(:, i)))
%     
%     %Max error of 10 folds
%     sprintf('Max error = %f', max(average_error31(:, i)))
%     
%     %Min error of 10 folds
%     sprintf('Min error = %f', min(average_error31(:, i)))
% end
% 
% % average_error1 =
% % 
% %     1.1939    1.8615    2.2245
% %     0.6065    0.4092    0.2973
% %     0.6034    0.4084    0.2973
% %     0.6168    0.4076    0.2999
% %     0.6051    0.4076    0.2995
% %     0.6096    0.4210    0.3011
% %     0.6143    0.4182    0.3002
% %     0.6017    0.4075    0.2940
% %     0.6009    0.4159    0.2986
% %     0.6056    0.4125    0.3050
% 
% % For k = 10 value
% % Average error = 0.665762
% % Max error = 1.193863
% % Min error = 0.600867
% % 
% % For k = 50 value
% % Average error = 0.556931
% % Max error = 1.861502
% % Min error = 0.407526
% % 
% % For k = 100 value
% % Average error = 0.491745
% % Max error = 2.224548
% % Min error = 0.293973
% 
% % -----------------------------------------------
% % Swapping R and W and using Regularized WNMFRule
% % -----------------------------------------------
% 
% 3 values of lambda
% 
% lambda = [0.01,0.1,1];
% 
% Training rating and weight matrices
% Swap values from previous part
% R_Actual32 = w31;
% w32 = R_Actual31;
% Matrix to store predicted values
% R_Predicted32= zeros(size(R_Actual32, 1), size(R_Actual32, 2), 3, 3); % [No. of Users, No. of movies, 3 values of k, 3 values of lambda]
% Matrix to store average error
% average_error32 = zeros(10, 3, 3); % [10 folds, 3 values of k, 3 values of lambda]
% 
% Training and Testing Data as per 10 - fold Cross Validation
% for j = 1:10
%     
%     Perform factorization, calculate error, and store predicted values
%     for n = 1:3 % For 3 values of lambda
%         for m = 1:3 % For 3 values of k
%             Perform factorization (using regularized WNMFRule) and get the Predicted values
%             [U3,V3,numIter,tElapsed,finalResidual]=mywnmfrulefor5(R_Actual32,k(m), w32, lambda(n));
% 
%             Predicted values
%             product1 = U3*V3;
% 
%             calculate the sum of absolute error and predicted matrix storage
%             sum_of_errors = 0;
%             for i=((j-1)*10000 + 1):j*10000
%                 Store indices as row3 and col3
%                 row3 = user_id(random_indices(i));
%                 col3 =  movie_id(random_indices(i));
% 
%                 Compute absolute error and add it to sum
%                 sum_of_errors = sum_of_errors + abs(product1(row3, col3) - R_Actual32(row3, col3));
% 
%                 Store the predicted value
%                 R_Predicted32(row3, col3, m, n) = product1(row3, col3);
%             end
% 
%             determine the average absolute error 
%             average_error32(j, m, n) = sum_of_errors / 10000
%         end
%     end
% end
% 
% Now print the average, maximum and minimum amongst 10 tests for each of the k values
% for j = 1:3 % For 3 values of lambda
%     for i = 1:3 % For 3 values of k
%         sprintf('For k = %d and lambda = %d', k(i), lambda(j))
% 
%         Average error of 10 folds
%         sprintf('Average error = %f', mean(average_error32(:, i, j)))
% 
%         Max error of 10 folds
%         sprintf('Max error = %f', max(average_error32(:, i, j)))
% 
%         Min error of 10 folds
%         sprintf('Min error = %f', min(average_error32(:, i, j)))
%     end
% end
% 
% % average_error1 =
% % 
% %     1.1939    1.8615    2.2245
% %     0.6065    0.4092    0.2973
% %     0.6034    0.4084    0.2973
% %     0.6168    0.4076    0.2999
% %     0.6051    0.4076    0.2995
% %     0.6096    0.4210    0.3011
% %     0.6143    0.4182    0.3002
% %     0.6017    0.4075    0.2940
% %     0.6009    0.4159    0.2986
% %     0.6056    0.4125    0.3050
% 
% % For k = 10 value
% % Average error = 0.665762
% % Max error = 1.193863
% % Min error = 0.600867
% % 
% % For k = 50 value
% % Average error = 0.556931
% % Max error = 1.861502
% % Min error = 0.407526
% % 
% % For k = 100 value
% % Average error = 0.491745
% % Max error = 2.224548
% % Min error = 0.293973
% 
% % -----------------------------------------------
% % Part 4
% % -----------------------------------------------
% 
% % Define Threshold value
% threshold = linspace(0.1, 5.0, 50); % Create a column vector with different threshold values
% 
% % Create precision and recall matrix 
% precision_total41 = zeros(size(threshold, 2), 3); % [No. of Threshold Values, 3 different values of k]
% recall_total41 = zeros(size(threshold, 2), 3); % [No. of Threshold Values, 3 different values of k]
% 
% % Precision and Recall computation
% for i=1:size(threshold, 2) % Loop over various threshold values
%     for m = 1:3 % For k values = 10, 50, 100
%         precision_total41(i, m)=length(find((R_Predicted31(:, :, m)>threshold(1,i)) & (R>3)))/length(find(R_Predicted31(:, :, m)>threshold(1,i)));
%         recall_total41(i, m)=length(find((R_Predicted31(:, :, m)>threshold(1,i)) & (R>3)))/length(find(R>3));
%     end
% end
% 
% Plot precision over recall values
% figure;
% plot(recall_total41(:, 1), precision_total41(:, 1), 'r', recall_total41(:, 2), precision_total41(:, 2), 'g', recall_total41(:, 3), precision_total41(:, 3), 'b')
% title('Precision versus Recall (Unregularized wnmf)')
% xlabel('Recall')
% ylabel('Precision')
% legend('k = 10', 'k = 50', 'k = 100')
% 
% % Plot precision over threshold values
% figure;
% plot(threshold(:), precision_total41(:, 1), 'r', threshold(:), precision_total41(:, 2), 'g', threshold(:), precision_total41(:, 3), 'b')
% title('Precision versus Threshold (Unregularized wnmf)')
% xlabel('Threshold')
% ylabel('Precision')
% legend('k = 10', 'k = 50', 'k = 100')
% 
% % Plot recall over threshold values
% figure;
% plot(threshold(:), recall_total41(:, 1), 'r', threshold(:), recall_total41(:, 2), 'g', threshold(:), recall_total41(:, 3), 'b')
% title('Recall versus Threshold (Unregularized wnmf)')
% xlabel('Threshold')
% ylabel('Recall')
% legend('k = 10', 'k = 50', 'k = 100')
% 
% % precision_total1 =
% % 
% %     0.5539    0.5541    0.5542
% %     0.5539    0.5543    0.5546
% %     0.5540    0.5546    0.5552
% %     0.5541    0.5549    0.5555
% %     0.5542    0.5552    0.5559
% %     0.5544    0.5556    0.5563
% %     0.5546    0.5559    0.5568
% %     0.5550    0.5564    0.5572
% %     0.5555    0.5571    0.5581
% %     0.5563    0.5589    0.5603
% %     0.5582    0.5642    0.5680
% %     0.5593    0.5664    0.5710
% %     0.5603    0.5690    0.5735
% %     0.5617    0.5714    0.5762
% %     0.5634    0.5737    0.5784
% %     0.5654    0.5765    0.5808
% %     0.5675    0.5793    0.5834
% %     0.5703    0.5823    0.5862
% %     0.5735    0.5864    0.5897
% %     0.5775    0.5919    0.5953
% %     0.5830    0.6019    0.6108
% %     0.5886    0.6095    0.6200
% %     0.5954    0.6177    0.6286
% %     0.6026    0.6266    0.6374
% %     0.6116    0.6366    0.6460
% %     0.6219    0.6476    0.6554
% %     0.6327    0.6602    0.6665
% %     0.6456    0.6752    0.6800
% %     0.6603    0.6932    0.6976
% %     0.6775    0.7160    0.7239
% %     0.6978    0.7548    0.7813
% %     0.7190    0.7857    0.8224
% %     0.7408    0.8161    0.8597
% %     0.7625    0.8451    0.8911
% %     0.7850    0.8710    0.9171
% %     0.8065    0.8955    0.9390
% %     0.8291    0.9169    0.9559
% %     0.8506    0.9351    0.9694
% %     0.8719    0.9496    0.9792
% %     0.8903    0.9628    0.9862
% %     0.9067    0.9712    0.9905
% %     0.9218    0.9793    0.9941
% %     0.9352    0.9854    0.9959
% %     0.9484    0.9896    0.9972
% %     0.9612    0.9930    0.9981
% %     0.9688    0.9956    0.9990
% %     0.9780    0.9970    0.9994
% %     0.9838    0.9989    0.9997
% %     0.9856    0.9993    0.9997
% %     0.9889    0.9996    0.9999
% 
% % recall_total1 =
% % 
% %     0.9999    0.9997    0.9995
% %     0.9999    0.9989    0.9984
% %     0.9997    0.9983    0.9973
% %     0.9996    0.9975    0.9954
% %     0.9995    0.9967    0.9930
% %     0.9993    0.9957    0.9905
% %     0.9990    0.9943    0.9876
% %     0.9987    0.9927    0.9840
% %     0.9984    0.9908    0.9799
% %     0.9980    0.9884    0.9753
% %     0.9973    0.9850    0.9707
% %     0.9966    0.9822    0.9653
% %     0.9957    0.9790    0.9599
% %     0.9949    0.9751    0.9546
% %     0.9938    0.9709    0.9487
% %     0.9926    0.9668    0.9430
% %     0.9910    0.9621    0.9376
% %     0.9892    0.9571    0.9323
% %     0.9873    0.9520    0.9276
% %     0.9848    0.9478    0.9231
% %     0.9823    0.9424    0.9191
% %     0.9795    0.9372    0.9156
% %     0.9759    0.9324    0.9124
% %     0.9711    0.9279    0.9096
% %     0.9656    0.9234    0.9075
% %     0.9589    0.9193    0.9054
% %     0.9502    0.9147    0.9035
% %     0.9401    0.9103    0.9020
% %     0.9281    0.9045    0.9002
% %     0.9130    0.8977    0.8978
% %     0.8946    0.8894    0.8949
% %     0.8726    0.8787    0.8902
% %     0.8447    0.8645    0.8831
% %     0.8118    0.8455    0.8737
% %     0.7730    0.8212    0.8590
% %     0.7251    0.7905    0.8378
% %     0.6747    0.7522    0.8063
% %     0.6185    0.7069    0.7657
% %     0.5560    0.6534    0.7116
% %     0.4927    0.5894    0.6419
% %     0.4245    0.5111    0.5379
% %     0.3582    0.4458    0.4632
% %     0.2963    0.3878    0.4071
% %     0.2406    0.3344    0.3618
% %     0.1908    0.2857    0.3221
% %     0.1480    0.2408    0.2850
% %     0.1117    0.2001    0.2479
% %     0.0831    0.1631    0.2114
% %     0.0592    0.1287    0.1717
% %     0.0402    0.0977    0.1313
% 
% % -----------------------------------------------
% % Swapping R and W and using Regularized WNMFRule
% % -----------------------------------------------
% 
% % Define Threshold value CHANGE THIS
threshold = linspace(0.01, 1, 50); % Create a column vector with different threshold values

% Create precision and recall matrix 
precision_total42 = zeros(size(threshold, 2), 3, 3); % [No. of Threshold Values, 3 different values of k, 3 different values of lambda]
recall_total42 = zeros(size(threshold, 2), 3, 3); % [No. of Threshold Values, 3 different values of k, 3 different values of lambda]

% Precision and Recall computation
for i=1:size(threshold, 2) % Loop over various threshold values
    for n = 1:3 % For lambda values = 0.01, 0.1, 1
        for m = 1:3 % For k values = 10, 50, 100
            precision_total42(i, m, n)=length(find((R_Predicted32(:, :, m, n)>threshold(1, i)) & (R>3)))/length(find(R_Predicted32(:, :, m, n)>threshold(1,i)));
            recall_total42(i, m, n)=length(find((R_Predicted32(:, :, m, n)>threshold(1, i)) & (R>3)))/length(find(R>3));
        end
    end
end

for n = 1:3 % For 3 values of lambda
    % Plot precision over recall values
    figure;
    plot(recall_total42(:, 1), precision_total42(:, 1), 'r', recall_total42(:, 2), precision_total42(:, 2), 'g', recall_total42(:, 3), precision_total42(:, 3), 'b')
    str = sprintf('Precision versus Recall (Regularized wnmf) for lambda = %d', lambda(n));
    title(str)
    xlabel('Recall')
    ylabel('Precision')
    legend('k = 10', 'k = 50', 'k = 100')

    % Plot precision over threshold values
    figure;
    plot(precision_total42(:, 1), threshold(:), 'r', precision_total42(:, 2), threshold(:), 'g', precision_total42(:, 3), threshold(:), 'b')
    str = sprintf('Precision versus Threshold (Regularized wnmf) for lambda = %d', lambda(n));
    title(str)
    xlabel('Threshold')
    ylabel('Precision')
    legend('k = 10', 'k = 50', 'k = 100')

    % Plot recall over threshold values
    figure;
    plot(recall_total42(:, 1), threshold(:), 'r', recall_total42(:, 2), threshold(:), 'g', recall_total42(:, 3), threshold(:), 'b')
    str = sprintf('Recall versus Threshold (Regularized wnmf) for lambda = %d', lambda(n));
    title(str)
    xlabel('Threshold')
    ylabel('Recall')
    legend('k = 10', 'k = 50', 'k = 100')
end
% 
% %--------------------------------------------------------------------------
% % Part 5
% %--------------------------------------------------------------------------
% % 3 values of lambda
% lambda = [0.01,0.1,1];
% 
% % Build matrix R and w (weights matrix)
% R5 = R;
% w5 = w;
% 
% % 3 values of k
% k = [10, 50, 100];
% 
% least_squared_error51 = zeros(3,3);
% 
% for i=1:3
%     [U51,V51,numIter,tElapsed,finalResidual511]=mywnmfrulefor5(R5, k(1), w5, lambda(i));
%     [U52,V52,numIter,tElapsed,finalResidual512]=mywnmfrulefor5(R5, k(2), w5, lambda(i));
%     [U53,V53,numIter,tElapsed,finalResidual513]=mywnmfrulefor5(R5, k(3), w5, lambda(i));
%     
%     % Compute and Store the errors
%     least_squared_error51(1,i) = finalResidual511; 
%     least_squared_error51(2,i) = finalResidual512; 
%     least_squared_error51(3,i) = finalResidual513; 
% end
% 
% % least_squared_error51 =
% % 
% %   245.6486  247.3199  247.7667
% %   174.7564  174.1914  179.2799
% %   132.4688  132.9680  138.3598
% 
% %-----------------------------------------------
% % Swapping R and W
% %-----------------------------------------------
% R52 = w5;
% w52 = R5;
% 
% least_squared_error52 = zeros(3,3);
% 
% for i=1:3
%     [U51,V51,numIter,tElapsed,finalResidual521]=mywnmfrulefor5(R52, k(1), w52, lambda(i));
%     [U52,V52,numIter,tElapsed,finalResidual522]=mywnmfrulefor5(R52, k(2), w52, lambda(i));
%     [U53,V53,numIter,tElapsed,finalResidual523]=mywnmfrulefor5(R52, k(3), w52, lambda(i));
%     
%     % Compute and Store the errors
%     least_squared_error52(1,i) = finalResidual521; 
%     least_squared_error52(2,i) = finalResidual522; 
%     least_squared_error52(3,i) = finalResidual523; 
% end
% 
% % least_squared_error52 =
% % 
% %    17.2589   16.2644   22.9115
% %    26.5806   25.7656   27.8253
% %    22.9265   23.9091   26.1655
% 
% %--------------------------------------------------------------------------
% % Part 6
% %--------------------------------------------------------------------------
% % STEP 1:
% % Convert R to a 0 - 1 Matrix, and use ratings as weights
% R6 = zeros(943, 1682);
% R_Predicted6 = R_Predicted32;
% 
% % STEP 2:
% % Find the top L movies
% L = 20;
% lambda = [0.01, 0.1, 1]; % All Lambda values
% % Create a matrix to store the top L movie IDs for each user, for k = 10,
% % 50, 100
% recommended_movies = zeros(size(R6, 1), L, 3, 3); %[No. of Users, L, k=[10, 50, 100], lamda = [0.01, 0.1, 1]]
% 
% % Variable for precision
% individual_precision = zeros(943, 3, 3); % [943 users, 3 values of k, 3 values of lambda]
% % Variable for Hit Rate
% individual_hit_rate = zeros(943, 20, 3, 3); % [943 users, 20 values of L, 3 values of k, 3 values of lambda]
% % Variable for False Alarm Rate
% individual_false_alarm_rate = zeros(943, 20, 3, 3); % [943 users, 20 values of L, 3 values of k, 3 values of lambda]
% threshold = 0.4;
% 
% % Find the top L movies
% for n = 1:3 % For 3 values of lambda
%     for m = 1:3 % For 3 different k values = 10, 50, 100
%         for i=1:size(R6, 1) % For each user
%             % Sort the Predicted Ratings for each user
%             [sorted_ratings, sorted_ratings_indices] = sort(R_Predicted6(i, :, m, n), 'descend');
%             
%             for L = 1:1:20 % Now pick the top L movies, for which we KNOW the ratings
%                 recommended_movies = zeros(size(R6, 1), L, 3, 3); % [No. of Users, L, k=[10, 50, 100], lambda = [0.01, 0.1, 1]]
%                 recommended_movies_counter = 1; % Initialize counter to 1
%                 for j=1:size(sorted_ratings_indices, 2)
%                     % Check if the ith user has given a rating for the jth movie
%                     if R2(i, sorted_ratings_indices(j)) == 1
%                         % add into the list of recommended movies
%                         recommended_movies(i, recommended_movies_counter, m, n) = (sorted_ratings_indices(j));
%                         recommended_movies_counter = recommended_movies_counter + 1;
%                     end
%                     
%                     % Check if we have already found the top L movies
%                     if recommended_movies_counter == L+1
% 
%                         % Calculate TP, TN, FP and FN
%                         true_positive = length(find((R_Predicted6(i, recommended_movies(i, :, m, n), m, n)>threshold) & (R(i, recommended_movies(i, :, m, n))>2)));
%                         true_negative = length(find((R_Predicted6(i, recommended_movies(i, :, m, n), m, n)<=threshold) & (R(i, recommended_movies(i, :, m, n))<=2)));
%                         false_positive = length(find((R_Predicted6(i, recommended_movies(i, :, m, n), m, n)>threshold) & (R(i, recommended_movies(i, :, m, n))<=2)));
%                         false_negative = length(find((R_Predicted6(i, recommended_movies(i, :, m, n), m, n)<=threshold) & (R(i, recommended_movies(i, :, m, n))>2)));
%                         
%                         % Calculate individual precision for the ith user, and
%                         % mth value of k and nth value of lambda ONLY IF L
%                         % = 5
%                         if L == 5
%                             
%                             individual_precision(i, m, n)= true_positive / length(find((R_Predicted6(i, recommended_movies(i, :, m, n), m, n)>threshold)));
%                         end
%                         
%                         % Calculate Hit Rate
%                         individual_hit_rate(i, L, m, n) = true_positive / (true_positive + false_negative);
%                         
%                         if true_positive == 0 & false_negative == 0
%                             individual_hit_rate(i, L, m, n) = 0;
%                         end
%                         
%                         % Calculate False Alarm Rate
%                         individual_false_alarm_rate(i, L, m, n) = false_positive / (false_positive + true_negative);
%                         
%                         if false_positive == 0 & true_negative == 0
%                             individual_false_alarm_rate(i, L, m, n) = 0;
%                         end
%                         
%                         break;
%                     end
%                 end 
%             end
%         end
%     end
% end
% 
% mean_hit_rate = zeros(20, 3, 3); % [20 values of L, 3 values of k, 3 values of lambda]
% mean_false_alarm_rate = zeros(20, 3, 3); % [20 values of L, 3 values of k, 3 values of lambda]
% 
% % Display average precision values and store average hit rate and false alarm rate for each value of k and lambda
% for n = 1:3 % For 3 values of lambda
%     for m = 1:3 % For 3 different k values = 10, 50, 100
%         mean(individual_precision(:, m, n))
%         for L = 1:1:20
%             mean_hit_rate(L, m, n) = mean(individual_hit_rate(:, L, m, n));
%             mean_false_alarm_rate(L, m, n) = mean(individual_false_alarm_rate(:, L, m, n));
%         end
%     end
% end
% 
% % Mean Precision values:
% %     0.7934
% %     0.7665
% %     0.7779
% %     0.7994
% %     0.7981
% %     0.8019
% %     0.8335
% %     0.8335
% %     0.8208
% 
% % Plot Mean Hit Rate and False Alarm Rate
% for n = 1:3 % For 3 values of lambda
%     figure;
%     plot(mean_false_alarm_rate(:, 1, n), mean_hit_rate(:, 1, n), 'r', mean_false_alarm_rate(:, 2, n), mean_hit_rate(:, 2, n), 'g', mean_false_alarm_rate(:, 3, n), mean_hit_rate(:, 3, n), 'b')
%     str = sprintf('Hit Rate versus False Alarm Rate (Regularized wnmf) for Lambda = %d', lambda(n));
%     title(str)
%     xlabel('False Alarm Rate')
%     ylabel('Hit Rate')
%     legend('k = 10', 'k = 50', 'k = 100')
% end

% mean_hit_rate(:,:,1) =
% 
%     0.7508    0.6957    0.7073
%     0.9162    0.9003    0.8940
%     0.9682    0.9682    0.9565
%     0.9809    0.9841    0.9809
%     0.9894    0.9936    0.9905
% 
% 
% mean_hit_rate(:,:,2) =
% 
%     0.7752    0.7253    0.7169
%     0.9374    0.9088    0.9236
%     0.9745    0.9618    0.9661
%     0.9915    0.9862    0.9883
%     0.9947    0.9958    0.9968
% 
% 
% mean_hit_rate(:,:,3) =
% 
%     0.8155    0.8006    0.7869
%     0.9608    0.9449    0.9332
%     0.9841    0.9852    0.9767
%     0.9936    0.9936    0.9905
%     0.9979    0.9947    0.9958

% mean_false_alarm_rate(:,:,1) =
% 
%     0.2492    0.3043    0.2927
%     0.3849    0.4464    0.4401
%     0.4719    0.5504    0.5154
%     0.5408    0.6204    0.5854
%     0.6013    0.6670    0.6299
% 
% 
% mean_false_alarm_rate(:,:,2) =
% 
%     0.2248    0.2747    0.2831
%     0.3669    0.4030    0.4019
%     0.4613    0.4804    0.4952
%     0.5366    0.5472    0.5514
%     0.6023    0.5970    0.5970
% 
% 
% mean_false_alarm_rate(:,:,3) =
% 
%     0.1845    0.1994    0.2131
%     0.3160    0.3351    0.3415
%     0.4040    0.4263    0.4327
%     0.4740    0.4793    0.4910
%     0.5398    0.5217    0.5483
