"""Run in the actual OmicVerse environment before installing scissor_py."""
import importlib.util
import sys
from offline_guard import install
blocked = install()
assert importlib.util.find_spec('scissor_py') is None
import omicverse as ov
cls = ov.single.Scissor
assert 'scissor_py' not in sys.modules
assert 'Scissor' in ov.single.__all__
model = cls(None, bulk=None, phenotype=None, family='binomial', graph_key='network')
try:
    model.run()
except ImportError as error:
    assert 'supplied scissor-py wheel' in str(error)
else:
    raise AssertionError('Missing backend was not reported')
assert model.result_ is None
print({'passed':True, 'omicverse':ov.__version__, 'path':ov.__file__,
       'missing_backend':True, 'blocked_network_calls':len(blocked)})
