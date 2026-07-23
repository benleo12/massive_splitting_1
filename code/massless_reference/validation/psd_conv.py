import pySecDec as psd
from pySecDec import LoopIntegralFromPropagators
li = LoopIntegralFromPropagators(
  propagators=['l**2 - 1','(l+k1)**2 - 0','(l+k1+k2)**2 - 1'], loop_momenta=['l'],
  external_momenta=['k1','k2'],
  replacement_rules=[('k1*k1','1'),('k2*k2','-1/5'),('k1*k2','(-3-1-(-1/5))/2')],
  regulators=['eps'], dimensionality='4-2*eps')
print("U =", li.U)
print("F =", li.F)
print("Gamma/prefactor =", li.Gamma_factor)
print("exponentF =", li.exponent_F, " exponentU =", li.exponent_U)
print("measure/extra =", getattr(li,'preliminary_F',None))
