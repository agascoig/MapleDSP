
# Holton DSP, p.456 Example 7.11 : Inverse DFT frequency-sampled filter
# Type-3 FIR bandpass filter
restart;
with(SignalProcessing):
N  := 17:
wc := Pi/2: # center frequency
dw := 0.5*Pi:
Nh := iquo(N-1,2):
w  := Vector(Nh, i -> 2*Pi*(i-1)/N): # frequency vector
Ak := Vector(Nh, i -> `if`(is(abs(w[i]-wc)<=0.5*dw), 1, 0)): # amplitude samples
Hk := Vector(Nh, i -> I*Ak[i]*exp(-I*Pi*(i-1)*(N-1)/N)): # complex spectrum
Hk_mirror := Vector(Nh, i -> conjugate(Hk[Nh - i + 1])):# mirror for linear phase FIR
H := Concatenate(1, Hk, Hk_mirror): # impulse response
h := InverseFFT(Hk, normalization=full): # frequency response
HF := FilterFrequencyResponse(h, size=1024, output=record, fftnormalization=none):
HF[magnitudeplot];
HF[phaseplot];
freq := HF[frequencies];
response := HF[response];
plot(freq, abs(response), title = "Absolute response magnitude", color = blue, gridlines = true, labels = ["Normalized Frequency (&times; &pi; rad/sample)", typeset("|H(", omega, ")|")], labeldirections = [horizontal, vertical]);
NULL;
