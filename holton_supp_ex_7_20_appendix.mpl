                            "maple in ~/.mapleinit"

with(LinearAlgebra): 
assume(n,integer): 
#&coloneq;0.45&#960;:
#&coloneq;0.55&#960;:
D2 := ptmp;
W := ptmp;
d := ptmp;
d := "1/(2 Pi) (&int;)[-omega[p]]^(omega[p])(e)^(omega n &ImaginaryI;) "" &DifferentialD;omega";
d2 := convert(d,trig);
d3 := subs(sin(omega__p*n)/Pi/n = omega__p/Pi*sinc(omega__p*n),d2);
w := `/`(1,2*Pi)*int(W(omega)*exp(I*omega*n),omega = -Pi .. Pi): 
w := omega__p/Pi*sinc(omega__p*n)+K*(delta(n)-omega__s/Pi*sinc(omega__s*n));
sinc := ptmp;

