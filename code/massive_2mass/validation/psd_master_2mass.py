#!/usr/bin/env python3
# pySecDec reference values for the 15 masters of the TWO-mass squark amplitude.
# Usage: python3 psd_master_2mass.py <master_id> <s> <t> <u> <m1sq> <m2sq>
import sys, subprocess, shutil
import pySecDec as psd
from pySecDec import LoopIntegralFromPropagators, loop_package
from fractions import Fraction
def fmt(x):
    fr=Fraction(x).limit_denominator(10**6)
    return str(fr.numerator) if fr.denominator==1 else '(%d/%d)'%(fr.numerator,fr.denominator)

mid = sys.argv[1]
s,t,u,m1sq,m2sq = [float(x) for x in sys.argv[2:7]]
Q = s+t+u-m1sq-m2sq                                   # Q^2 = s+t+u-m1^2-m2^2
V = {'0':0.0, 'm1':m1sq, 'm2':m2sq, 'sM':s,'tM':t,'uM':u,'QM':Q}

def bub(psq, ma, mb):
    return (['l**2 - (%s)'%fmt(V[ma]), '(l+p1)**2 - (%s)'%fmt(V[mb])], ['l'], ['p1'],
            [('p1*p1', V[psq])])
def tri(p1,p2,p3, ma,mb,mc):
    P1,P2,P3 = V[p1],V[p2],V[p3]
    return (['l**2 - (%s)'%fmt(V[ma]), '(l+k1)**2 - (%s)'%fmt(V[mb]), '(l+k1+k2)**2 - (%s)'%fmt(V[mc])],
            ['l'], ['k1','k2'], [('k1*k1',P1), ('k2*k2',P2), ('k1*k2',(P3-P1-P2)/2)])
def box(p1,p2,p3,p4, s12,s23, ma,mb,mc,md):
    P1,P2,P3,P4,S12,S23 = V[p1],V[p2],V[p3],V[p4],V[s12],V[s23]
    return (['l**2 - (%s)'%fmt(V[ma]), '(l+k1)**2 - (%s)'%fmt(V[mb]),
             '(l+k1+k2)**2 - (%s)'%fmt(V[mc]), '(l+k1+k2+k3)**2 - (%s)'%fmt(V[md])],
            ['l'], ['k1','k2','k3'],
            [('k1*k1',P1),('k2*k2',P2),('k3*k3',P3),
             ('k1*k2',(S12-P1-P2)/2),('k2*k3',(S23-P2-P3)/2),('k1*k3',(P4-S12-S23+P2)/2)])

M = {
 # bubbles
 'B0_m1_0_m1': bub('m1','0','m1'),   'B0_m2_0_m2': bub('m2','0','m2'),
 'B0_s_m1_m2': bub('sM','m1','m2'),  'B0_Q_m1_m2': bub('QM','m1','m2'),
 'B0_t_0_m1':  bub('tM','0','m1'),   'B0_u_0_m2':  bub('uM','0','m2'),
 # triangles
 'C0_0_m1_t':  tri('0','m1','tM',  'm1','m1','0'),
 'C0_0_m2_u':  tri('0','m2','uM',  'm2','m2','0'),
 'C0_0_s_Q_a': tri('0','sM','QM',  'm1','m1','m2'),
 'C0_0_s_Q_b': tri('0','sM','QM',  'm2','m2','m1'),
 'C0_m1_m2_s': tri('m1','m2','sM', 'm1','0','m2'),
 'C0_m1_u_Q':  tri('m1','uM','QM', 'm1','0','m2'),
 'C0_m2_t_Q':  tri('m2','tM','QM', 'm2','0','m1'),
 # boxes
 'D0_a': box('0','uM','m1','sM', 'm2','QM', 'm2','m2','0','m1'),
 'D0_b': box('m1','tM','QM','sM','0','m2',  'm1','0','m1','m2'),
}
prop, loop, ext, rr = M[mid]
li = LoopIntegralFromPropagators(propagators=prop, loop_momenta=loop,
     external_momenta=ext, replacement_rules=[(a,fmt(b)) for a,b in rr],
     regulators=['eps'], dimensionality='4-2*eps')
name='m_'+mid
shutil.rmtree(name, ignore_errors=True)
loop_package(name=name, loop_integral=li, requested_orders=[3])
subprocess.run(['make','-C',name,'-j4'],check=True,stdout=subprocess.DEVNULL)
from pySecDec.integral_interface import IntegralLibrary
lib=IntegralLibrary(name+'/'+name+'_pylink.so'); lib.use_Qmc(verbosity=0)
_,_,wp=lib()
print("PSDRESULT %s s=%g t=%g u=%g m1sq=%g m2sq=%g Q2=%g :: %s"%(mid,s,t,u,m1sq,m2sq,Q,wp))
