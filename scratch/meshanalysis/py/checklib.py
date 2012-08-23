
import operator
import inspect

_failures = []
def check_failure(msg):
    global _failures
    _failures.append(msg)
    print msg

def check_failures():
    global _failures
    return _failures

def check_op(name, op, a, b):
    res = op(a, b)
    if not res:
        #frame = inspect.currentframe()
        stack = inspect.stack()
        code = None
        for k,l in enumerate(stack):
            code = l[4][0].strip()
            if code.startswith('check'):
                break
        msg = "Check failed: %s %s %s, from code:\n    %s" % (a, name, b, code)
        check_failure(msg)
    return res

def check(a):
    return check_op("", lambda x,y: bool(x), a, "")

def check_eq(a, b):
    return check_op("==", operator.eq, a, b)

def check_lt(a, b):
    return check_op("<", operator.lt, a, b)

def check_gt(a, b):
    return check_op(">", operator.gt, a, b)

def check_le(a, b):
    return check_op("<=", operator.le, a, b)

def check_ge(a, b):
    return check_op(">=", operator.ge, a, b)
