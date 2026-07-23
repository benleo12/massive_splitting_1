# C0(0,s,Q^2; mq2,mq2,mq2) at s=-1, Q2=-5, mq2=1.  k1^2=0,k2^2=s=-1,k1.k2=(Q2-s)/2=-2
import subprocess
import pySecDec as psd
from pySecDec import LoopIntegralFromPropagators, loop_package
li=LoopIntegralFromPropagators(
  propagators=['l**2-1','(l+k1)**2-1','(l+k1+k2)**2-1'], loop_momenta=['l'],
  external_momenta=['k1','k2'],
  replacement_rules=[('k1*k1','0'),('k2*k2','-1'),('k1*k2','-2')],
  regulators=['eps'], dimensionality='4-2*eps')
loop_package(name='c0sq',loop_integral=li,requested_orders=[3])
subprocess.run(['make','-C','c0sq','-j4'],check=True,stdout=subprocess.DEVNULL)
from pySecDec.integral_interface import IntegralLibrary
lib=IntegralLibrary('c0sq/c0sq_pylink.so');lib.use_Qmc(verbosity=0)
_,_,wp=lib();print("PSD_C0(0,s=-1,Q2=-5;mq2=1):",wp)
