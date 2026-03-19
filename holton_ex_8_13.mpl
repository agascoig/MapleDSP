
# Holton p. 560 Example 8.13 Type-II discrete-time inverse Chebyshev lowpass filter
# Design a discrete-time inverse (Type-II) Chebyshev lowpass filter using the bilinear-transformation method using the same parameters as the previous example: ωp = 0.5π, δp = -3 dB, ωs = 0.7π and δs = -40 dB.
restart;
with(ArrayTools);
with(PolynomialTools);
with(SignalProcessing);
assume(z, real);
omega__p := 0.50000*Pi;
omega__s := 0.70000*Pi;
d__p := -3;
d__s := -40;
WpT := 2*tan(omega__p/2);
WsT := 2*tan(omega__s/2);
e__s := evalf(sqrt(10^(-d__s/10) - 1));
e__p := evalf(sqrt(10^(-d__p/10) - 1));
x := arccosh(e__s/e__p);
N := ceil(x/arccosh(WsT/WpT));
WcT := WsT;
NULL;
NULL;
k := Vector(N, i -> i - 1);
psi := `+`~(2.00000*k, 1)*Pi/(2*N);
phi := arcsinh(e__s)/N;
pkT := `/`~(WcT, -sin~(psi)*sinh(phi) + cos~(psi)*cosh(phi)*I);
psi := Remove(psi, floor(N/2) + 1);
zkT := I . WcT /~ cos~(psi);
pls := `+`~(2, pkT) /~ `-`~(2, pkT);
zrs := Concatenate(1, `+`~(2, zkT) /~ `-`~(2, zkT), -1);
b := expand(Re~(mul((zkT[i] - 2)/zkT[i], i = 1 .. NumElems(zkT))*mul(pkT[i]/(pkT[i] - 2), i = 1 .. NumElems(pkT)))*Re~(mul(z - zr, zr in zrs)));
a := Re~(expand(mul(z - p, p in pls)));
b2 := Vector([seq(coeff(b, z, i), i = degree(b, z) .. 0, -1)]);
a2 := Vector([seq(coeff(a, z, i), i = degree(a, z) .. 0, -1)]);
HF := FilterFrequencyResponse([b2], [a2], size = 1024, output = record, fftnormalization = none);
HF[magnitudeplot];
HF[phaseplot];
freq := HF[frequencies];
response := HF[response];
plot(freq, abs(response), title = "Absolute response magnitude", color = blue, gridlines = true, labels = ["Normalized Frequency (&times; &pi; rad/sample)", typeset("|H(", omega, ")|")], labeldirections = [horizontal, vertical]);
