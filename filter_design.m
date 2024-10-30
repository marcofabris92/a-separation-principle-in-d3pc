% This AUXILIARY function is used to design a low pass filter, whose
% transfer function coefficients are stored in HLP
% Given
% - Fst: Stop-Band frequency
% - Ts: sampling time

% Invoked by:
% - data_generation.m, to filter the input signal u_long if requested
% Invkoes: none


function HLP = filter_design(Ts,Fst)

%Fst = Stop-Band frequency  1/pi
%Fp = Pass band frequency

Fp = Fst - 0.2/(2*pi);

% tolerances
% in pass band (deltap)
% and stop band (delats)

deltap = 0.01; 
deltas = 0.001;

% Tolerances expressed in dB

Ap = 20*log10(1+deltap)-20*log10(1-deltap); % (approximate)
Ast = -20*log10(deltas);



SPEC = 'Fp,Fst,Ap,Ast';
LP = fdesign.lowpass(SPEC,Fp,Fst,Ap,Ast,1/Ts);

HLP = design(LP,'equiripple');




% plot frequency response: the red line is the ideal low pass filter with
% tolerances and transition band depicted; the blue line is the realized 
% frequency response

% fvtool(HLP)

end