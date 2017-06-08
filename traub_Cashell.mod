TITLE Calcium dynamics Traub

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (mS) = (millisiemens)
}

NEURON {
    THREADSAFE
    : note - every variable accessible in NEURON will be having the suffix _CaShell
    
    SUFFIX CaShell
    USEION ca READ ica WRITE cai
    RANGE phi
    : GLOBAL tnmax, tlmax
}

PARAMETER {
    phi = 33.2223   : value at soma
    beta_x = 0.075
    cai0 = 0
}

STATE {
    cai
}

ASSIGNED {
    ica
}

BREAKPOINT {
    SOLVE states METHOD cnexp
}

INITIAL {
   cai = cai0
}

DERIVATIVE states {  
    cai' = - phi*ica - beta_x*cai
}

