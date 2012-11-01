
import sys
from math import round
from decimal import Decimal

import numpy as np

def msg(msg='', *args):
    print >>sys.stderr, msg%args if args else msg
    sys.stderr.flush()

class FauxRepr(object):
    def __init__(self, op):
        self.op = op
    def __repr__(self):
        return str(self.op)

def float_fix(op, places=2):
    factor = 10 ** places
    return FauxRepr(Decimal(int(round(op * factor))) / factor)

def epsilon_eq(a, b, epsilon=None):
    if np.isnan(a) or np.isnan(b):
        return False
    if epsilon is None:
        epsilon = 1e-3 * (abs(a) + abs(b))
    assert epsilon >= 0, 'bad epsilon %r' % (epsilon,)
    return abs(a-b) <= epsilon

def all_epsilon_eq(seq, val, epsilon=None):
    seq = np.asarray(seq)
    if np.isnan(seq).any() or np.isnan(val):
        return False
    if epsilon is None:
        epsilon = 1e-3 * (abs(seq.mean() + abs(val)))
    assert epsilon >= 0, 'bad epsilon %r' % (epsilon,)
    return (abs(seq-val) <= epsilon).all()
