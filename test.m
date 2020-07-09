clear all
syms m S_theta I_theta K_h K_theta b p q v PI S e
Ms = [m         S_theta;
      S_theta   I_theta];

Ks = [K_h    0;
       0     K_theta];

   
Ca = (q/v) * [-S*2*PI      -S*b*2*PI;
              S*b*PI           0    ];
          
Ka = [0     -2*PI*S;
      0     e*2*b*S*2*PI];
  
res = Ms*p^2 - Ca*p + Ks - q * Ka;
det_ = det(res);
[t, coeff] = coeffs(det_, p);
t(1)
t(2)
t(3)
t(4)
t(5)