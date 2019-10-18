
function ins_list = inertial_navigation_system(accel_readings, gyro_readings, pt_0, vt_0, e_angles)
    
    f_e = accel_readings;
    omega_iee = gyro_readings;
%     time = time_stamp;
    
    gravity = [0; 0; 9.81]; % m/s^2
    
    pos0 = pt_0;
    vel0 = vt_0;
    euler_angles_0 = e_angles;
    % euler_angles_0 = [yaw_0; ptc_0; rol_0];

    dt = 0.01;

    % x_t = (position, velocity, orientation)
    x0 = [pos0; vel0; euler_angles_0];
    x_t = x0;

    mu_lat = 47.6*pi/180; % deg
    mu_long = -122.3*pi/180; % deg

    % yaw -> pitch -> roll
    % psi -> theta -> phi
    % euler_angles = (yaw, pitch, roll)

    % omega_icc = omege_itc
    % earths rotation
    omega_icc = [0; 0; (2*pi / 86400)]; % rad/s
    omega_itc = omega_icc;
    
    ins_list = zeros(numel(10000, 9));
    
    for m1 = 1:10000
        
        ins_list(1, m1) = x_t(1);
        ins_list(2, m1) = x_t(2);
        ins_list(3, m1) = x_t(3);
        ins_list(4, m1) = x_t(4);
        ins_list(5, m1) = x_t(5);
        ins_list(6, m1) = x_t(6);
        ins_list(7, m1) = x_t(7);
        ins_list(8, m1) = x_t(8);
        ins_list(9, m1) = x_t(9);

        f_e_t = f_e(:, m1);
        omega_iee_t = omega_iee(:, m1);
        x_t = rk4_step(x_t, f_e_t, omega_iee_t);

    end

    function x_dot = mech_eq_rhs(x_t, f, o)

        vel_t = x_t(4:6);
        psi_m = x_t(7);
        theta_m = x_t(8);
        phi_m = x_t(9);

        Rct = [cos(-(0.5*pi + mu_lat)) 0 -sin(-(0.5*pi + mu_lat));
                             0         1            0            ;
               sin(-(0.5*pi + mu_lat)) 0 cos(-(0.5*pi + mu_lat))] * ...
               [cos(mu_long) sin(mu_long) 0;
                -sin(mu_long) cos(mu_long) 0;
                0               0           1];

        H321e = [-sin(theta_m) 0 1;
                 sin(phi_m)*cos(theta_m) cos(phi_m) 0;
                 cos(phi_m)*cos(theta_m) -sin(phi_m) 0];

        Rte_f = [1 0 0; 0 cos(phi_m) sin(phi_m); 0 -sin(phi_m) cos(phi_m);];
        Rte_s = [cos(theta_m) 0 -sin(theta_m); 0 1 0; sin(theta_m) 0 cos(theta_m)];
        Rte_t = [cos(psi_m) sin(psi_m) 0; -sin(psi_m) cos(psi_m) 0; 0 0 1];
             
        Rte = Rte_f*Rte_s*Rte_t;
        
        Ret = Rte';

        omega_itt = Rct*omega_icc;


        pos_dot = vel_t;
        vel_dot = Ret*f + gravity - 2*cross(omega_itt, vel_t);
        e_angle_dot = H321e\(o - Rte*Rct*omega_itc); 

        x_dot = [pos_dot; vel_dot; e_angle_dot];

    end

    function x_t_plus_dt = rk4_step(x_t, f, o)

        k1 = dt * mech_eq_rhs(x_t, f, o);
        k2 = dt * mech_eq_rhs(x_t + 0.5*k1, f, o);
        k3 = dt * mech_eq_rhs(x_t + 0.5*k2, f, o);
        k4 = dt * mech_eq_rhs(x_t + k3, f, o);

        x_t_plus_dt = x_t + (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4;

    end

end



