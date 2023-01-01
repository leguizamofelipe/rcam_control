function [xdot] = rcam_model(x, u)
    %------------------------ STATE AND CONTROL ---------------------------
    % Extract state vector
    x1 = x(1); % u
    x2 = x(2); % v
    x3 = x(3); % w
    x4 = x(4); % p
    x5 = x(5); % q
    x6 = x(6); % r
    x7 = x(7); % phi
    x8 = x(8); % theta
    x9 = x(9); % psi
    
    % Extract control inputs
    u1 = u(1); % d_A: Aileron
    u2 = u(2); % d_T: Stabilizer
    u3 = u(3); % d_r: Rudder
    u4 = u(4); % d_th1: Throttle 1
    u5 = u(5); % d_th2: Throttle 2
    
    %-------------------------- CONSTANTS ---------------------------------
    m = 120000; % AC mass (kg)

    cbar = 6.6; % Mean aerodynamic chord (m)
    lt = 24.8;  % Distance by AC of tail and body (m)
    S = 260;    % Wing platform area (m^2)
    St = 64;    % Tail platform area (m^2)
    
    % x, y, z position of CG (m)
    Xcg = 0.23*cbar; 
    Ycg = 0;
    Zcg = 0.1*cbar;
    
    % x, y, z position of aerodynamic center (m)
    Xac = 0.12*cbar;
    Yac = 0;
    Zac = 0;

    % x, y, z position of engine 1 force (m)
    Xapt1 = 0;
    Yapt1 = -7.94;
    Zapt1 = -1.9;

    % x, y, z position of engine 2 force (m)
    Xapt2 = 0;
    Yapt2 = 7.94;
    Zapt2 = -1.9;

    % Other constants
    rho = 1.225;                    % Air density (kg/m^3)
    g = 9.81;                       % Gravitational acceleration (m/s^2)
    deps_da = 0.25;                 % Change in downwash wrt to alpha (rad/rad)
    alpha_L0 = -11.5*pi/180;        % Zero lift angle of attack (rad)
    n = 5.5;                        % Slope of linear region of lift slope
    a3 = -768.5;                    % Coefficient of alpha^3
    a2 = 609.2;                     % Coefficient of alpha^2
    a1 = -155.2;                    % Coefficient of alpha^1
    a0 = 15.212;                    % Coefficient of alpha^0
    alpha_switch = 14.5*(pi/180);   % alpha where slope is no longer linear (rad)

    %------------------- CONTROL LIMITS AND SATURATION --------------------
    % Implement in Simulink

    %------------------- INTERMEDIATE VARIABLES ---------------------------
    % Calculate airspeed
    Va = sqrt(x1^2 + x2^2 + x3^2);
    
    % Calculate alpha and beta
    alpha = atan2(x3, x1);
    beta = asin(x2/Va);

    % Calculate dynamic pressure
    Q = 0.5*rho*Va^2;

    % Define vectors wbe_b and V_b
    wbe_b = [x4; x5; x6];
    V_b = [x1; x2; x3];

    %------------------ AERODYNAMIC FORCE COEFFICIENTS --------------------
    % Calculate the Cl of the wing body
    if alpha <= alpha_switch
        % Linear region
        Cl_wb = n*(alpha-alpha_L0);
    else
        Cl_wb = a3*alpha^3 + a2*alpha^2 + a1*alpha + a0;
    end

    % Calculate Cl_t
    epsilon = deps_da*(alpha-alpha_L0);
    alpha_t = alpha - epsilon + u2 + 1.3*x5*lt/Va;
    Cl_t = 3.1*(St/S)*alpha_t;

    % Total lift force 
    Cl = Cl_wb + Cl_t;
    
    % Total drag force (neglecting tail)
    Cd = 0.13 + 0.07*(5.5*alpha + 0.654)^2;

    % Calculate sideforce
    Cy = -1.6*beta + 0.24*u3;
    
    % Calculate the actual dimensional forces. These are in F_s (stability axis)
    FA_s = [-Cd*Q*S;
             Cy*Q*S;
            -Cl*Q*S];

    % Rotate these forces to F_b (body axis)
    C_bs = [cos(alpha) 0 -sin(alpha);
            0 1 0;
            sin(alpha) 0 cos(alpha)];

    FA_b = C_bs*FA_s;

    %--------------- AERODYNAMIC MOMENT COEFFICIENT ABOUT AC --------------
    % Calculate the moments in Fb. Define eta, dCmdx, and dCMdu
    eta11 = -1.4*beta;
    eta21 = -0.59 - (3.1*(St*lt)/(S*cbar))*(alpha-epsilon);
    eta31 = (1-alpha*(180/(15*pi)))*beta;

    eta = [eta11;
           eta21;
           eta31;];

    dCm_dx = (cbar/Va)*[-11 0 5;
                        0 (-4.03*(St*lt^2)/(S*cbar^2)) 0;
                        1.7 0 -11.5];

    dCm_du = [-0.6 0 0.22;
              0 (-3.1*(St*lt)/(S*cbar)) 0;
              0 0 -0.63];

    Cm_ac_b = eta + dCm_dx*wbe_b + dCm_du*[u1; u2; u3];

    MA_ac_b = Cm_ac_b*Q*S*cbar;

    %----------------- AERODYNAMIC MOMENT ABOUT CG ------------------------
    rcg_b = [Xcg; Ycg; Zcg];
    rac_b = [Xac; Yac; Zac];
    MA_cg_b = MA_ac_b + cross(FA_b, rcg_b-rac_b);

    %------------------- ENGINE FORCE AND MOMENT --------------------------
    % Thrust from each engine
    F1 = u4*m*g;
    F2 = u5*m*g;

    FE1_b = [F1;0;0];
    FE2_b = [F2;0;0];
    
    FE_b = FE1_b + FE2_b;

    % Now engine moment due to offset of engine thrust from CG
    mew1 = [Xcg - Xapt1;
            Yapt1 - Ycg;
            Zcg - Zapt1];

    mew2 = [Xcg - Xapt2;
            Yapt2 - Ycg;
            Zcg - Zapt2];

    MEcg1_b = cross(mew1, FE1_b);
    MEcg2_b = cross(mew2, FE2_b);

    MEcg_b = MEcg1_b + MEcg2_b;

    %------------------- GRAVITY EFFECTS ----------------------------------
    % Calculate gravitational forces in the body frame. No moment about CG
    g_b = [  -g*sin(x8);
            g*cos(x8)*sin(x7);
            g*cos(x8)*cos(x7)];

    Fg_b = m*g_b;

    %------------------------ INERTIA MATRIX ------------------------------
    Ib = m*[40.07 0 -2.0923;
    0 64 0;
    -2.0923 0 99.92];

    invIb = (1/m) * [0.0249836 0 0.00523151;
    0 0.015625 0;
    0.000523151 0 0.010019];

    F_b = Fg_b + FE_b + FA_b;
    xlto3dot = (1/m) * F_b - cross(wbe_b,V_b);
    
    Mcg_b = MA_cg_b + MEcg_b;
    x4to6dot = invIb * (Mcg_b - cross(wbe_b,Ib * wbe_b));
    
    H_phi = [   1 sin(x7)*tan(x8) cos(x7)*tan(x8)
                0 cos(x7) -sin(x7);
                0 sin(x7)/cos(x8) cos(x7)/cos(x8)];
            
    x7to9dot = H_phi * wbe_b;
    
    xdot = [xlto3dot
            x4to6dot
            x7to9dot];

end