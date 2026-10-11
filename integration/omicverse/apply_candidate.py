"""Apply this local candidate to the exact pinned OmicVerse source checkout.

This never clones, fetches, installs, commits, or publishes. It refuses to
replace unrelated/modified files and changes only single/__init__.py plus the
new single/_scissor.py module. Run against an explicit private source copy.
"""
import argparse
import hashlib
from pathlib import Path

COMMIT = 'f9e1db8875d750d7f679514243a86fd2aa37d770'
# Filled from the untouched official pinned source, checked before mutation.
INIT_SHA256 = '64475abd116f77a5836607cf5cd1a50f8ea887b787cf74dbe05ab96b4844cd7a'


def apply(root):
    original = root/'omicverse/single/__init__.py'
    module = root/'omicverse/single/_scissor.py'
    candidate = Path(__file__).parent/'omicverse/single/_scissor.py'
    text = original.read_text()
    changed = text.replace('from ._cnv import CNV\n', 'from ._cnv import CNV\nfrom ._scissor import Scissor\n', 1)
    changed = changed.replace('__all__ = [\n', "__all__ = [\n    'Scissor',\n", 1)
    if hashlib.sha256(original.read_bytes()).hexdigest() != INIT_SHA256:
        expected = text.replace('from ._scissor import Scissor\n', '', 1).replace("    'Scissor',\n", '', 1)
        if hashlib.sha256(expected.encode()).hexdigest() == INIT_SHA256 and module.exists() and module.read_bytes() == candidate.read_bytes():
            print('Exact candidate already applied; no changes.')
            return
        raise SystemExit('Refusing modified/unrecognized single/__init__.py; expected pinned '+COMMIT)
    if module.exists():
        raise SystemExit('Refusing to overwrite an existing single/_scissor.py')
    if changed == text or changed.count('from ._scissor import Scissor') != 1:
        raise SystemExit('Registration anchor mismatch')
    # Both checks precede the two small local writes. Restore original on error.
    module.write_bytes(candidate.read_bytes())
    try: original.write_text(changed)
    except BaseException:
        module.unlink()
        original.write_text(text)
        raise
    print('Applied local Scissor candidate to '+str(root))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkout', type=Path)
    apply(parser.parse_args().checkout)
