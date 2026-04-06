                            "maple in ~/.mapleinit"

# Holton DSP, p.456 Example 7.11 : Inverse DFT frequency-sampled filter
# Type-3 FIR bandpass filter
restart: 
N := 17: 
eps := .1e-3: 
wc := `/`(Pi,2): # center frequency;
dw := `*`(.5,Pi): 
Nh := iquo(N-1,2): 
w := Vector(Nh+1,i -> 2*Pi*(i-1)/N): # frequency vector;
Ak := Vector(Nh+1,i -> `if`(is(abs(w[i]-wc) <= .5*dw),1,eps)): # amplitude samples;
Hk := Vector(Nh+1,i -> I*Ak[i]*exp(-I*Pi*(i-1)*(N-1)/N)): # complex spectrum for Type-3 FIR;
Hk_mirror := Vector(Nh,i -> conjugate(Hk[Nh+2-i])): # mirror for linear phase FIR;
H := Concatenate(1,Hk,Hk_mirror): # impulse response;
stemplot(abs(H));
h := Re(InverseFFT(H,normalization = full)): # frequency response;
stemplot(abs(h));
HF := FilterFrequencyResponse(h,size = 1024,output = record,fftnormalization = none): 
HF[magnitudeplot];
HF[phaseplot];
freq := HF[frequencies]: 
response := HF[response]: 
plot(freq,abs(response),title = "Absolute response magnitude",color = blue,gridlines = true,labels = ["Normalized Frequency (&times; &pi; rad/sample)", typeset("|H(",omega,")|")],labeldirections = [horizontal, vertical]);

