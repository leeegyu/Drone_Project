function X_DOT = Dynamics_GyL_Ver02(STATE, CONTROL, Prop)

MASS = Prop.MASS;
g = Prop.g;
J = Prop.J;
Cn2b = Prop.Cn2b;
l = Prop.l;
k_F = Prop.k_F;
k_M = Prop.k_M;

P = STATE(1:3, 1);
V = STATE(4:6, 1);
AngVel = STATE(7:9, 1);
Q = STATE(10:13, 1);

p = AngVel(1);
q = AngVel(2);
r = AngVel(3);

w = CONTROL;

Qmat = [0 -p -q -r;
    p 0 r -q;
    q -r 0 p;
    r q -p 0];

Rmat = [-k_F,             -k_F,             -k_F,            -k_F;
         k_F * sqrt(2)/2, -k_F * sqrt(2)/2, -k_F * sqrt(2)/2, k_F * sqrt(2)/2; 
        -k_F * sqrt(2)/2, -k_F * sqrt(2)/2,  k_F * sqrt(2)/2, k_F * sqrt(2)/2;
         k_M,              k_M,             -k_M,            -k_M];

FMmatrix = Rmat * (w.^2);

Fz = FMmatrix(1);
L = FMmatrix(2);
M = FMmatrix(3);
N = FMmatrix(4);

F = [0; 0; Fz];

%6-DoF eq.
Pos_dot = (Cn2b')*V;
Vel_dot = (1/MASS)*F - cross(AngVel,V) + Cn2b*g;
AngVel_dot = J\(M - cross(AngVel, J*AngVel));
Q_dot = 0.5*Qmat*Q;

X_DOT = [Pos_dot; Vel_dot; AngVel_dot; Q_dot];

end

