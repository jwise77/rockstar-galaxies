#ifndef _UNIVERSAL_CONSTANTS_H_
#define _UNIVERSAL_CONSTANTS_H_

#define Gc (4.30117902e-9) //Actually, Gc * (Msun / Mpc) in (km/s)^2
#define CRITICAL_DENSITY 2.77519737e11 // 3H^2/8piG in (Msun / h) / (Mpc / h)^3
#define VMAX_CONST 6.55833746e-5 //sqrt(G*(Msun/h)/(Mpc/h)) in km/s
#define RMAX_TO_RS 2.1626        /* Rvmax / Rs */
#define RS_CONSTANT 0.216216595 /* ln(1+2.1626)/2.1626 - 1/(1+2.1626) */

#define HUBBLE_TIME_CONVERSION 9.77813952e9 /* (100 km/s/Mpc)^-1 to years */
#define ENERGY_PER_KELVIN_PER_UNIT_MASS 0.00712668983 /* 1.5 * (1 / (0.75*amu + 0.25*4 amu))*k to ((km/s)^2 / Kelvin) */ 

#endif /* _UNIVERSAL_CONSTANTS_H_ */
