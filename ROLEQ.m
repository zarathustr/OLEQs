% Recursive Optimal Linear Estimator of Quaternion (ROLEQ)
% Author: Jin Wu, Zebo Zhou, Jinling Wang and Hassen Fourati
% e-mail: jin_wu_uestc@hotmail.com, klinsmann.zhou@gmail.com


function q = ROLEQ(Gyroscope, dt, Db, Dr, weights, q_last, iter)

    len = length(Db(1, :));
    
    W = zeros(4, 4);
    
    wx = Gyroscope(1);      wy = Gyroscope(2);      wz = Gyroscope(3);
    
    for i = 1 : len
        
        W = W + weights(i) * 0.5 * (W_matrix(Db(:, i), Dr(:, i)) + eye(4));
        
    end
    
    Omg = [0, -wx, -wy, -wz;
           wx, 0, wz, -wy;
           wy, -wz, 0, wx;
           wz, wy, -wx, 0] * dt * 0.5 + eye(4);

    q =  W^iter * Omg * q_last;

    q = q ./ norm(q);
end