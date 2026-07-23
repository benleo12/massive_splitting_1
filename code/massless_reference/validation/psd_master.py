#!/usr/bin/env python3
# pySecDec reference values for the 13 massive masters.
# Usage: python3 psd_master.py <master_id> <s> <t> <u> <mq2>
import sys, subprocess
import pySecDec as psd
from pySecDec import LoopIntegralFromPropagators, loop_package
from fractions import Fraction
def fmt(x):
    fr=Fraction(x).limit_denominator(10**6)
    return str(fr.numerator) if fr.denominator==1 else '(%d/%d)'%(fr.numerator,fr.denominator)

mid = sys.argv[1]
s,t,u,mq2 = [float(x) for x in sys.argv[2:6]]
Q = s+t+u-2*mq2
V = {'0':0.0, 'mq2':mq2, 'sM':s, 'tM':t, 'uM':u, 'QM':Q}     # invariant tokens -> numbers

def bub(psq, m1, m2):
    return (['l**2 - (%s)'%fmt(V[m1]), '(l+p1)**2 - (%s)'%fmt(V[m2])], ['l'], ['p1'],
            [('p1*p1', V[psq])])

def tri(p1,p2,p3, m1,m2,m3):
    P1,P2,P3 = V[p1],V[p2],V[p3]
    return (['l**2 - (%s)'%fmt(V[m1]), '(l+k1)**2 - (%s)'%fmt(V[m2]), '(l+k1+k2)**2 - (%s)'%fmt(V[m3])],
            ['l'], ['k1','k2'],
            [('k1*k1',P1), ('k2*k2',P2), ('k1*k2',(P3-P1-P2)/2)])

def box(p1,p2,p3,p4, s12,s23, m1,m2,m3,m4):
    P1,P2,P3,P4,S12,S23 = V[p1],V[p2],V[p3],V[p4],V[s12],V[s23]
    return (['l**2 - (%s)'%fmt(V[m1]), '(l+k1)**2 - (%s)'%fmt(V[m2]),
             '(l+k1+k2)**2 - (%s)'%fmt(V[m3]), '(l+k1+k2+k3)**2 - (%s)'%fmt(V[m4])],
            ['l'], ['k1','k2','k3'],
            [('k1*k1',P1),('k2*k2',P2),('k3*k3',P3),
             ('k1*k2',(S12-P1-P2)/2),('k2*k3',(S23-P2-P3)/2),
             ('k1*k3',(P4-S12-S23+P2)/2)])

m='mq2'; Z='0'
M = {
 'B0_t_0_mq': bub('tM',Z,m), 'B0_u_0_mq': bub('uM',Z,m), 'B0_mq2_0_mq': bub('mq2',Z,m),
 'B0_s_mq_mq': bub('sM',m,m), 'B0_Q_mq_mq': bub('QM',m,m),
 'C0_0_mq_t': tri('0','mq2','tM', m,m,Z), 'C0_0_mq_u': tri('0','mq2','uM', m,m,Z),
 'C0_0_s_Q': tri('0','sM','QM', m,m,m), 'C0_mq_mq_s': tri('mq2','mq2','sM', m,Z,m),
 'C0_mq_t_Q': tri('mq2','tM','QM', m,Z,m), 'C0_mq_u_Q': tri('mq2','uM','QM', m,Z,m),
 'D0_a': box('0','uM','mq2','sM', 'mq2','QM', m,m,Z,m),
 'D0_b': box('mq2','tM','QM','sM', '0','mq2', m,Z,m,m),
}
prop, loop, ext, rr = M[mid]
li = LoopIntegralFromPropagators(propagators=prop, loop_momenta=loop,
     external_momenta=ext, replacement_rules=[(a,fmt(b)) for a,b in rr],
     regulators=['eps'], dimensionality='4-2*eps')
name='m_'+mid
import shutil; shutil.rmtree(name, ignore_errors=True)   # regen fresh: kinematics are baked in at codegen
loop_package(name=name, loop_integral=li, requested_orders=[3])
subprocess.run(['make','-C',name,'-j4'],check=True,stdout=subprocess.DEVNULL)
from pySecDec.integral_interface import IntegralLibrary
lib=IntegralLibrary(name+'/'+name+'_pylink.so'); lib.use_Qmc(verbosity=0)
_,_,wp=lib()
print("PSDRESULT %s s=%g t=%g u=%g mq2=%g Q2=%g :: %s"%(mid,s,t,u,mq2,Q,wp))
