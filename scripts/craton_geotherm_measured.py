"""Section 24.7: rerun of the craton geotherm with the measured xenolith values (Section 24.6).
Same conductive model as craton_geotherm.py (exponential heat production, constant k).
Forward: Moho temperature at q0 = 40-47 mW/m2 with heat production 0.3-0.6 uW/m3.
Inverse: surface heat flow q0 and mantle heat flow q_m needed to give the measured Moho temperature.
Measured (Rudnick, Apen and Ringwood 2024 abstract): Siberia < 400 C, A0 ~0.3; Tanzania 500-600 C, A0 ~0.5.
Assumed (memory, unsourced): T0 = 15 C, decay depth hr 8-12 km, k 2.5-3.0 W/m/K, Moho depth 35-42 km."""
import itertools
import numpy as np

T0 = 15.0
def moho_T(q0, a0, hr, k, z):
    q0 *= 1e-3; a0 *= 1e-6; hr *= 1e3; z *= 1e3
    return T0 + (q0 - a0 * hr) * z / k + a0 * hr**2 / k * (1 - np.exp(-z / hr))

def q0_for(T, a0, hr, k, z):            # linear in q0, so invert directly
    t0 = moho_T(0.0, a0, hr, k, z)
    t1 = moho_T(1.0, a0, hr, k, z)
    return (T - t0) / (t1 - t0)

HR, K, ZM = (8.0, 10.0, 12.0), (2.5, 3.0), (35.0, 39.0, 42.0)

print("FORWARD: Moho temperature (C) with heat production 0.3-0.6 uW/m3")
for q0 in (40.0, 44.0, 47.0):
    T = np.array([moho_T(q0, a, h, k, z) for a in (0.3, 0.45, 0.6) for h in HR for k in K for z in ZM])
    print(f"  q0 {q0:.0f} mW/m2: median {np.median(T):.0f}, range {T.min():.0f}-{T.max():.0f}, "
          f"share >= 600 C {100*(T>=600).mean():.0f}%")

def inverse(name, Tlo, Thi, a0):
    qs, qms = [], []
    for T, h, k, z in itertools.product((Tlo, Thi), HR, K, ZM):
        q = q0_for(T, a0, h, k, z); qs.append(q); qms.append(q - a0 * h)
    qs, qms = np.array(qs), np.array(qms)
    print(f"  {name:22s} Moho {Tlo:.0f}-{Thi:.0f} C, A0 {a0}: surface q0 {qs.min():4.1f}-{qs.max():4.1f} "
          f"(median {np.median(qs):4.1f}) mW/m2; mantle heat flow q_m {qms.min():4.1f}-{qms.max():4.1f} mW/m2")

print("\nINVERSE: heat flow needed to give the measured Moho temperature")
inverse("Siberia (upper bound)", 300, 400, 0.3)
inverse("Tanzania craton", 500, 600, 0.5)
print("  (Siberian measured surface heat flow 19 mW/m2 is called 'impossibly low' by the source.)")

print("\nPrediction check: Moho temperature if surface q0 were only 19 mW/m2 (Siberian measurement) with A0 0.3:")
T = np.array([moho_T(19.0, 0.3, h, k, z) for h in HR for k in K for z in ZM])
print(f"  median {np.median(T):.0f} C, range {T.min():.0f}-{T.max():.0f}")
