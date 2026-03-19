
# Holton DSP, p. 454 Example 7.10: Inverse DFT frequency-sampled filter Type-1 lowpass filter
restart;
with(SignalProcessing):
N  := 17:
wc := Pi/2:
Nh := iquo(N,2):
w := Vector(Nh, i -> 2*Pi*(i-1)/N): # frequency samples
Ak := Vector(Nh, i -> `if`(is(w[i] < wc), 1, 0)): # amplitude samples
Hk := Vector(Nh,
      i -> Ak[i]*exp(-I*Pi*(i-1)*(N-1)/N)): # complex spectrum
Hk_mirror := Vector(Nh, i -> conjugate(Hk[Nh - i + 1])): # mirror for linear phase FIR
H := Concatenate(1, Hk, Hk_mirror):
h := InverseFFT(Hk, normalization=full): # impulse response
HF := FilterFrequencyResponse(h, size=1024, output=record, fftnormalization=none): # frequency response
HF[magnitudeplot];
HF[phaseplot];
freq := HF[frequencies];
response := HF[response];
plot(freq, abs(response), title = "Absolute response magnitude", color = blue, gridlines = true, labels = ["Normalized Frequency (&times; &pi; rad/sample)", typeset("|H(", omega, ")|")], labeldirections = [horizontal, vertical]);
NULL;
