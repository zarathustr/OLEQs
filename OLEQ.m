% Optimal Linear Estimator of Quaternion (OLEQ)
% Author: Jin Wu, Zebo Zhou, Jinling Wang and Hassen Fourati
% e-mail: jin_wu_uestc@hotmail.com, klinsmann.zhou@gmail.com


function [q, times] = OLEQ(Db, Dr, weights)

    len = length(Db(1, :));
    
    W = zeros(4, 4);
    
    for i = 1 : len
        
        W = W + weights(i) * W_matrix(Db(:, i), Dr(:, i));
        
    end
    
    last_q = [1; 0; 0; 0];
    q = [0; 0; 0; 1];
    times = 0;
    G = 0.5 * (W + eye(4));
    G = G * G;
    G = G * G;
    
    while(norm(q - last_q) > 1e-8 && times < 50)
        last_q = q;
        q = G * last_q;
        q = q ./ norm(q);
        times = times + 1;
    end

    q = q ./ norm(q);
end