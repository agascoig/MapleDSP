                            "maple in ~/.mapleinit"

# Parks McClellan Remez Exchange Method
# Use the Parks-McClellan algorithm to design a linear-phase, symmetric (Type-1) lowpass filter of odd length N.  From Holton DSP Supplementary.
Entered April 4th, 2026.
restart; 
omega__p := `*`(.25,Pi): 
omega__s := `*`(.5,Pi): 
delta__ps := 2: 
N := 13: 
# The response of the filter is given by the equation A = omega -> sum(a[m] cos(m omega),m=0..M).;
initial_omegas := Vector([0, omega__p, `/`(omega__p,2), omega__s, 1.25*omega__s, 1.5*omega__s, 1.75*omega__s, `*`(1.0,Pi)]): 
pm_iter := proc(N, omega__p,omega__s,delta__ps,omegas)
	local cosomegas, M, W, DD, c, delta, AB, L, AL, EST_A,
	S, error_freqs, error_A, freqs, k, i, A, omega, w, pi;

M := iquo(N-1,2):

cosomegas := cos~(omegas);

# Weighting
W := proc(omega)
	return `if`(is(0 <= omega and omega <= omega__p), 1.0, delta__ps);
end proc;

# Desired response
DD := proc(omega)
	return `if`(is(0 <= omega and omega <= omega__p), 1.0, 0.0);
end proc;

c := proc(k)
	local i; 
	return mul(seq(`if`(is(i <> k), 1/(cosomegas[k + 1] - cosomegas[i + 1]),
	NULL), i = 0 .. M + 1));
end proc;

# Current maximum error # TBD: not correct value
delta := add(c(k)*DD(omegas[k + 1]), k = 0 .. M + 1)/add((-1)^k*c(k)/W(omegas[k + 1]), k = 0 .. M + 1);

# Use (7.62) to compute values of A(w) at the extremal frequencies.
AB := seq(DD(omegas[k + 1]) - (-1)^k*delta/W(omegas[k + 1]), k = 0 .. M + 1);

# Use Lagrangian polynomial interpolation to determine A(w) from the M+2 values of A(w__k).
L := proc(w, k)
	local j;
	return mul(seq(`if`(is(j <> k), (cos(w) - cosomegas[j + 1])/(cosomegas[k + 1] - cosomegas[j + 1]), NULL), j = 0 .. M));
end proc;
AL := proc(w)
	local k, pos;
	return `if`(member(w, omegas, 'pos'), AB[pos], add(AB[k + 1]*L(w, k), k = 0 ..M));
end proc;

# Compute estimates for the given frequencies
EST_A := AL~(omegas);

# Now use suggestion in Hinamoto "Digital Filter Design and Realization" to compute error values for 16*(N-1)..
S := 16*(N - 1);
pi := evalf(Pi);
error_freqs := Vector([seq(`if`(is(omega__p < omega and omega < omega__s), NULL, omega), omega = 0 .. pi, 1.00000/(S*pi))]);
error_A := Vector(NumElems(error_freqs),0.0);
i := 1;
for w in error_freqs do
	error_A[i]:=abs(W(w)*(DD(w)-AL(w)));
	i := i+1;
	end do:
	
# Now compute A function for entire 0..pi interval
freqs := Vector(2*S,k->(k-1)*pi/(2*S));
A := Vector(NumElems(freqs),0.0);
i := 1;
for w in freqs do
	A[i]:=AL(w);
	i := i+1;
	end do:	
return delta,EST_A,error_freqs,error_A,freqs,A;
end proc:
(delta, EST_A, error_freqs, error_A, f, A) := pm_iter(N,omega__p,omega__s,delta__ps,initial_omegas): 
plot(f,20*log10~(abs(A)));
EST_A;

