import subprocess
import pySecDec as psd
from pySecDec import LoopIntegralFromPropagators, loop_package
li = LoopIntegralFromPropagators(
    propagators=['l**2', '(l+p1)**2 - 1'], loop_momenta=['l'],
    external_momenta=['p1'], replacement_rules=[('p1*p1','-2')],
    regulators=['eps'], dimensionality='4-2*eps')
loop_package(name='b0m2', loop_integral=li, requested_orders=[3])
subprocess.run(['make','-C','b0m2','-j4'],check=True,stdout=subprocess.DEVNULL)
from pySecDec.integral_interface import IntegralLibrary
lib=IntegralLibrary('b0m2/b0m2_pylink.so'); lib.use_Qmc(verbosity=0)
_,_,wp = lib()
print("PSD_B0(p2=-2;0,1):", wp)
