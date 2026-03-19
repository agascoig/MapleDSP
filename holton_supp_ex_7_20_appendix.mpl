
with(LinearAlgebra);
assume(n, integer);
NULL;
NULL;
D2 := omega -> piecewise(abs(omega) < omega__p, 1, 0);
W := omega -> piecewise(abs(omega) < omega__p, 1, omega__p <= abs(omega) and abs(omega) < omega__s, 0, K);
d := n -> 1/2*Int(D2(omega)*W(omega)*exp^(omega*n*I), omega = -Pi .. Pi)/Pi;
d := 1/(2*Pi)*int(exp(omega*n*I), omega = -omega__p .. omega__p);
d2 := convert(d, trig);
d3 := subs(sin(omega__p*n)/(Pi*n) = omega__p*sinc(omega__p*n)/Pi, d2);
w := 1/(2*Pi)*int(W(omega)*exp(omega*n*I), omega = -Pi .. Pi);
w := omega__p*sinc(omega__p*n)/Pi + K*(delta(n) - omega__s*sinc(omega__s*n)/Pi);
sinc := x -> piecewise(x = 0, 1, sin(x)/x);
