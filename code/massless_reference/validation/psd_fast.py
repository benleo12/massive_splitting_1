#!/usr/bin/env python3
# FAST path (1-mass): compile each master ONCE with kinematics as pySecDec
# real_parameters, then evaluate any point by a sub-second library call.
#   build one:   python3 psd_fast.py build   <id>
#   eval one:    python3 psd_fast.py eval     <id> <s> <t> <u> <mq2>
#   eval all 13: python3 psd_fast.py evalall  <s> <t> <u> <mq2>
# real_parameters order: (s, t, u, mq2).  Q^2 = s+t+u-2*mq2.
import sys, shutil
V = {'0':'0','mq2':'mq2','sM':'s','tM':'t','uM':'u','QM':'(s+t+u-2*mq2)'}
RP = ['s','t','u','mq2']
def half(a,b,c): return '((%s)-(%s)-(%s))/2'%(V[a],V[b],V[c])
def bub(psq,ma,mb):
    return (['l**2 - (%s)'%V[ma], '(l+p1)**2 - (%s)'%V[mb]], ['l'], ['p1'], [('p1*p1', V[psq])])
def tri(p1,p2,p3, ma,mb,mc):
    return (['l**2 - (%s)'%V[ma], '(l+k1)**2 - (%s)'%V[mb], '(l+k1+k2)**2 - (%s)'%V[mc]],
            ['l'], ['k1','k2'], [('k1*k1',V[p1]), ('k2*k2',V[p2]), ('k1*k2', half(p3,p1,p2))])
def box(p1,p2,p3,p4, s12,s23, ma,mb,mc,md):
    return (['l**2 - (%s)'%V[ma], '(l+k1)**2 - (%s)'%V[mb],
             '(l+k1+k2)**2 - (%s)'%V[mc], '(l+k1+k2+k3)**2 - (%s)'%V[md]],
            ['l'], ['k1','k2','k3'],
            [('k1*k1',V[p1]),('k2*k2',V[p2]),('k3*k3',V[p3]),
             ('k1*k2',half(s12,p1,p2)),('k2*k3',half(s23,p2,p3)),
             ('k1*k3','((%s)-(%s)-(%s)+(%s))/2'%(V[p4],V[s12],V[s23],V[p2]))])
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
IDS = list(M.keys())
def build(mid):
    from pySecDec import LoopIntegralFromPropagators, loop_package
    import subprocess
    prop, loop, ext, rr = M[mid]
    li = LoopIntegralFromPropagators(propagators=prop, loop_momenta=loop, external_momenta=ext,
         replacement_rules=rr, regulators=['eps'], dimensionality='4-2*eps')
    name='m_'+mid; shutil.rmtree(name, ignore_errors=True)
    loop_package(name=name, loop_integral=li, requested_orders=[3], real_parameters=RP)
    subprocess.run(['make','-C',name,'-j4'],check=True,stdout=subprocess.DEVNULL)
    print("built %s"%mid)
def load(mid):
    from pySecDec.integral_interface import IntegralLibrary
    lib=IntegralLibrary('m_'+mid+'/m_'+mid+'_pylink.so'); lib.use_Qmc(verbosity=0); return lib
def emit(mid, lib, s,t,u,mq2):
    _,_,wp = lib(real_parameters=[s,t,u,mq2]); Q=s+t+u-2*mq2
    print("PSDRESULT %s s=%g t=%g u=%g mq2=%g Q2=%g :: %s"%(mid,s,t,u,mq2,Q,wp))
if __name__=='__main__':
    mode=sys.argv[1]
    if mode=='build': build(sys.argv[2])
    elif mode=='eval':
        mid=sys.argv[2]; s,t,u,mq2=[float(x) for x in sys.argv[3:7]]; emit(mid, load(mid), s,t,u,mq2)
    elif mode=='evalall':
        s,t,u,mq2=[float(x) for x in sys.argv[2:6]]
        for mid in IDS: emit(mid, load(mid), s,t,u,mq2)
