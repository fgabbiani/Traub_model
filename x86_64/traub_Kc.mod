TITLE KC channel Traub

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (mS) = (millisiemens)
}

NEURON {
    THREADSAFE
    : note - every variable accessible in NEURON will be having the suffix _gKc
    
    SUFFIX gKc
    USEION ca READ cai
    USEION k READ ek WRITE ik
    RANGE gmax, g, i
    GLOBAL c_inf, tau_c
}

PARAMETER {
    gmax=0.03 (mho/cm2)  : value at soma
}

STATE {
    c
}

ASSIGNED {
    v (mV)
    ek (mV)
    ik		(mA/cm2)
    c_inf
    tau_c	(ms)
    g		(S/cm2)
    i
    cai
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g = gmax*c*xi(cai)
    i = g*(v-ek)
    ik = i
}

INITIAL {
    rates(v)
    c = c_inf
}

FUNCTION xi(cai) {
    xi = cai/250
    if ( xi > 1 ) {
	xi = 1
    }
}

FUNCTION alpha_c(v(mV)) {
    if ( v <= 50 ) {
	alpha_c = exp( (v - 10)/11 - (v - 6.5)/27 )/18.975
    } else {
	alpha_c = 2*exp( -(v - 6.5)/27 )
    }
}

FUNCTION beta_c(v(mV)) {
    if ( v<= 50 ) {
	beta_c = 2*exp( -(v - 6.5)/27 ) - alpha_c(v)
    } else {
	beta_c = 0
    }
}

DERIVATIVE states {  
    rates(v)
    c' = (c_inf - c)/tau_c
}

PROCEDURE rates(v (mV)) { :callable from hoc
    LOCAL alpha, beta
    TABLE c_inf, tau_c  : DEPEND vhalfn, vhalfl, tlmax, tnmax
    FROM -20 TO 130 WITH 750
    
    alpha = alpha_c(v)
    beta = beta_c(v)
    tau_c = 1/(alpha + beta)     
    c_inf = alpha*tau_c
}

