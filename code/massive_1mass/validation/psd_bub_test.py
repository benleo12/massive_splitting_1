#!/usr/bin/env python3
# pySecDec pipeline test: massive bubble B0(psq; mq2, mq2) at a Euclidean point.
# Generates + compiles + evaluates -> eps-series (eps^-1 .. eps^3).
import sys, subprocess, os
import pySecDec as psd
from pySecDec import LoopIntegralFromPropagators, loop_package

NAME = "bub_mm"
li = LoopIntegralFromPropagators(
    propagators       = ['l**2 - mq2', '(l + p1)**2 - mq2'],
    loop_momenta      = ['l'],
    external_momenta  = ['p1'],
    replacement_rules = [('p1*p1', 'psq')],
    regulators        = ['eps'],
    dimensionality    = '4-2*eps',
)
# requested_orders=[3]: expand to eps^3 beyond the leading pole
loop_package(name=NAME, loop_integral=li, requested_orders=[3],
             real_parameters=['psq', 'mq2'])
print("generated; compiling...", flush=True)
subprocess.run(['make', '-C', NAME, '-j4'], check=True,
               stdout=subprocess.DEVNULL)
print("compiled; evaluating...", flush=True)
from pySecDec.integral_interface import IntegralLibrary
lib = IntegralLibrary(f'{NAME}/{NAME}_pylink.so')
lib.use_Qmc(verbosity=0)
# Euclidean point: psq = -1, mq2 = 1
pref, no_pref, with_pref = lib(real_parameters=[-1.0, 1.0])
print("=== prefactor ===");        print(pref)
print("=== integral (with prefactor) ==="); print(with_pref)
