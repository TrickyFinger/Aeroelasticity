syms m S_theta I_theta K_h K_theta b p q v PI S
Ms = [m         S_theta;
      S_theta   I_theta];

Ks = [K_h    0;
       0     K_theta];

   
Ca = (q/v) * [-S*2*PI      -S*b*2*PI;
              S*b*PI           0    ];
          
Ka = [0     -2*PI*S;
      0     0.5*b*S*2*PI];
  
res = Ms*p^2 - Ca*p + Ks - Ka;
det = det(res);
[t, coeffs] = coeffs(det, p)