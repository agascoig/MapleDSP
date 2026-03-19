
# Holton DSP, p.454 Example 7.12 : Linear system frequency-sampled filter
# Type-1 FIR lowpass filter of Ex 7.10 using simultanous equations
restart;
with(SignalProcessing): with(LinearAlgebra): with(ArrayTools):
N  := 17:
wc := Pi/2: # center frequency
Nh := iquo(N-1,2):
k := Vector(Nh+1, i -> i-1)^%T: # row vector
w  := (2*Pi*k^%H)/N: # frequency samples, column vector
Ak := Vector(Nh+1, i -> `if`(is(w[i] < wc), 1, 0)): # amplitude samples, column vector
C := cos~(evalf~(w . k)):
a := LinearSolve(C, Ak): # outer product
h := Concatenate(1, 0.5*Reverse(a)[1..NumElems(a)-1], a[1], 0.5*a[2..]): # impulse response
H := FilterFrequencyResponse(h, size=1024, output=record, fftnormalization=none):
H[magnitudeplot];
H[phaseplot];
freq := H[frequencies];
response := H[response];
plot(freq, abs(response), title = "Absolute response magnitude", color = blue, gridlines = true, labels = ["Normalized Frequency (&times; &pi; rad/sample)", typeset("|H(", omega, ")|")], labeldirections = [horizontal, vertical]);
NULL;
