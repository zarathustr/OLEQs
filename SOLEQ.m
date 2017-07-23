% Sub-Optimal Linear Estimator of Quaternion (SOLEQ)
% Author: Jin Wu, Zebo Zhou, Jinling Wang and Hassen Fourati
% e-mail: jin_wu_uestc@hotmail.com, klinsmann.zhou@gmail.com


function q = SOLEQ(Db, Dr)

    n = length(Db(1,:));

    if n > 1
    
        V1 = V_matrix(Db(:,1), Dr(:,1));
        Vn = V_matrix(Db(:,n), Dr(:,n));
    
        P = eye(4);
        
        if n > 2
        
            for i = 2 : n - 1
                P = P * (W_matrix(Db(:,i), Dr(:,i)) + eye(4)) * 0.5;
            end
        
        end
        
        N=[
            0, 0, 0, 0;
            0, 0, 0, 0;
            V1(1,3), V1(2,3), V1(3,3), V1(4,3);
            V1(1,4), V1(2,4), V1(3,4), V1(4,4);
        ] ...
        * P * ...
        [
            0, 0, Vn(1,3), Vn(1,4);
            0, 0, Vn(2,3), Vn(2,4);
            0, 0, Vn(3,3), Vn(3,4);
            0, 0, Vn(4,3), Vn(4,4);
        ];
        
        n1 = N(3, 3);
        n2 = N(3, 4);
        n3 = N(4, 3);
        n4 = N(4, 4);

        u1 = n1 * n1 + n2 * n2;
        u2 = n1 * n3 + n2 * n4;        
        u4 = n3 * n3 + n4 * n4;

        sqrt_delta = sqrt(u1 * u1 - 2 * u1 * u4 + 4 * u2 * u2 + u4 * u4);
    
        g2 = 0.5 * (u1 - u4 + sqrt_delta) / u2;
    
        q = [V1(1,4) + V1(1,3) * g2;
             V1(2,4) + V1(2,3) * g2;
             V1(3,4) + V1(3,3) * g2;
             V1(4,4) + V1(4,3) * g2];
        
    else
        
        q = ((W_matrix(Db(:,1), Dr(:,1)) + eye(4)) * 0.5) * [1;0;0;0];
        
    end
    
    q = q./norm(q);
end