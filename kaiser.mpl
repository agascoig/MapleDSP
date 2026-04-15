Holton p. 446 Kaiser Lowpass Filter
restart;
                    "maple in ~/.mapleinit"

omega__c := 0.50000*Pi;
delta__omega := 0.10000*Pi;
delta__p := 0.01000;
delta__s := 0.01000;
omega__p := omega__c - delta__omega/2;
omega__s := omega__c + omega__omega/2;
delta := min(delta__p, delta__s);
A := -20*log10(delta);
beta := piecewise(50 < A, 0.11020*(A - 8.70000), 21 <= A and A <= 50, 0.58420*(A - 21)^0.40000 + 0.07886*(A - 21), 0);
f__c := omega__p/Pi;
M := ceil((A - 7.95000)/(2.28500*delta__omega));
M := M + (M mod 2);
n := Vector(M + 1, i -> i - 1/2*M - 1);
w := Vector(NumElems(n), i -> MTM:-besselj(0, beta *~ sqrt(1 - `*`~(2, (i - 1/2*M - 1)/M) ^~ 2))/MTM:-besselj(0, beta));
w := evalf~(w);
h := `*`~(w, f__c) *~ sinc~(`*`~(f__c, n));
stemplot(h);

HF := FilterFrequencyResponse(h, size = 1024, output = record, fftnormalization = none);
HF[magnitudeplot];

HF[phaseplot];
