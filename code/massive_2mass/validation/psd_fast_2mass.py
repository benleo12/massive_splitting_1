#!/usr/bin/env python3
# FAST path: compile each master ONCE with kinematics as pySecDec real_parameters,
# then evaluate any point by a sub-second library call (no recompile).
#
#   build one:   python3 psd_fast_2mass.py build   <id>
#   eval one:    python3 psd_fast_2mass.py eval     <id> <s> <t> <u> <m1sq> <m2sq>
#   eval all 15: python3 psd_fast_2mass.py evalall  <s> <t> <u> <m1sq> <m2sq>
#
# real_parameters order is ALWAYS  (s, t, u, m1sq, m2sq).  Q^2 = s+t+u-m1sq-m2sq.
import sys, shutil

# kinematics kept SYMBOLIC (strings in the real_parameters)
V = {'0':'0','m1':'m1sq','m2':'m2sq','sM':'s','tM':'t','uM':'u','QM':'(s+t+u-m1sq-m2sq)'}
RP = ['s','t','u','m1sq','m2sq']
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

M = {
 'B0_m1_0_m1': bub('m1','0','m1'),   'B0_m2_0_m2': bub('m2','0','m2'),
 'B0_s_m1_m2': bub('sM','m1','m2'),  'B0_Q_m1_m2': bub('QM','m1','m2'),
 'B0_t_0_m1':  bub('tM','0','m1'),   'B0_u_0_m2':  bub('uM','0','m2'),
 'C0_0_m1_t':  tri('0','m1','tM',  'm1','m1','0'),
 'C0_0_m2_u':  tri('0','m2','uM',  'm2','m2','0'),
 'C0_0_s_Q_a': tri('0','sM','QM',  'm1','m1','m2'),
 'C0_0_s_Q_b': tri('0','sM','QM',  'm2','m2','m1'),
 'C0_m1_m2_s': tri('m1','m2','sM', 'm1','0','m2'),
 'C0_m1_u_Q':  tri('m1','uM','QM', 'm1','0','m2'),
 'C0_m2_t_Q':  tri('m2','tM','QM', 'm2','0','m1'),
 'D0_a': box('0','uM','m1','sM', 'm2','QM', 'm2','m2','0','m1'),
 'D0_b': box('m1','tM','QM','sM','0','m2',  'm1','0','m1','m2'),
}
IDS = list(M.keys())

def build(mid):
    from pySecDec import LoopIntegralFromPropagators, loop_package
    import subprocess
    prop, loop, ext, rr = M[mid]
    li = LoopIntegralFromPropagators(propagators=prop, loop_momenta=loop, external_momenta=ext,
         replacement_rules=rr, regulators=['eps'], dimensionality='4-2*eps')
    name='m_'+mid
    shutil.rmtree(name, ignore_errors=True)
    loop_package(name=name, loop_integral=li, requested_orders=[3], real_parameters=RP)
    subprocess.run(['make','-C',name,'-j4'],check=True,stdout=subprocess.DEVNULL)
    print("built %s"%mid)

def load(mid):
    from pySecDec.integral_interface import IntegralLibrary
    lib=IntegralLibrary('m_'+mid+'/m_'+mid+'_pylink.so'); lib.use_Qmc(verbosity=0)
    # Guard: baked-in-kinematics libs expect 0 real_parameters; calling with 5 raises
    # inside a worker thread and the main thread DEADLOCKS.  Fail loudly instead.
    n = int(lib.info['number_of_real_parameters'])
    if n != len(RP):
        raise SystemExit(
            "ERROR: m_%s expects %d real_parameters but this fast path supplies %d (%s).\n"
            "       Rebuild it:  python3 psd_fast_2mass.py build %s" % (mid, n, len(RP), ','.join(RP), mid))
    return lib

def emit(mid, lib, s,t,u,m1sq,m2sq):
    _,_,wp = lib(real_parameters=[s,t,u,m1sq,m2sq])
    Q=s+t+u-m1sq-m2sq
    print("PSDRESULT %s s=%g t=%g u=%g m1sq=%g m2sq=%g Q2=%g :: %s"%(mid,s,t,u,m1sq,m2sq,Q,wp))

if __name__=='__main__':
    mode=sys.argv[1]
    if mode=='build':
        build(sys.argv[2])
    elif mode=='eval':
        mid=sys.argv[2]; s,t,u,m1sq,m2sq=[float(x) for x in sys.argv[3:8]]
        emit(mid, load(mid), s,t,u,m1sq,m2sq)
    elif mode=='evalall':
        s,t,u,m1sq,m2sq=[float(x) for x in sys.argv[2:7]]
        for mid in IDS: emit(mid, load(mid), s,t,u,m1sq,m2sq)
