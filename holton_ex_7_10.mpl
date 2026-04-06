                            "maple in ~/.mapleinit"

# Holton DSP, p. 454 Example 7.10: Inverse DFT frequency-sampled filter Type-1 lowpass filter
restart: 
N := 17: 
wc := `/`(Pi,2): 
Nh := iquo(N,2): 
eps := .1e-3: 
w := Vector(Nh+1,i -> 2*Pi*(i-1)/N): # frequency samples;
Ak := Vector(Nh+1,i -> `if`(is(w[i] < wc or 2*Pi-w[i] < wc),1,eps)): # amplitude samples;
Hk := Vector(Nh+1,i -> Ak[i]*exp(-I*Pi*(i-1)*(N-1)/N)): # complex spectrum;
Hk_mirror := Vector(Nh,i -> conjugate(Hk[Nh+2-i])): # mirror for linear phase FIR;
H := Concatenate(1,Hk,Hk_mirror): 
stemplot(abs(H));
h := Re(InverseFFT(H,normalization = full)): # impulse response;
HF := FilterFrequencyResponse(h,size = 1024,output = record,fftnormalization = none): # noncausal frequency response;
HF[magnitudeplot];
HF[phaseplot];
freq := HF[frequencies]: 
response := HF[response]: 
plot(freq,abs(response),title = "Absolute response magnitude",color = blue,gridlines = true,labels = ["Normalized Frequency (&times; &pi; rad/sample)", typeset("|H(",omega,")|")],labeldirections = [horizontal, vertical]);

