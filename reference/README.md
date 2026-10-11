# Golden fixtures and development provenance

This directory contains saved synthetic expected outputs, SHA-256 records and
upstream license/source notices. It is not a runtime backend. The Python package
never imports or executes anything here.

- `fixtures/`: all 22 original R-generated JSON fixtures, required by R-free CI.
- `review/`: two additional synthetic random-path/CV sweep JSON files.
- `manifest.json`: unchanged original 31-file manifest, including the hashes of
  nine upstream source/notice files. Six R/C++ source files now live only in the
  archive; their hashes are retained, not presented as locally rechecked files.
- `PROVENANCE.json`: exact archive and upstream commits.
- `upstream/`: original GPL license, DESCRIPTION and README attribution notices.
- `sessionInfo.txt`: historical oracle environment, not runtime requirements.

## Full development archive

The complete original audit, all 28 R/C++ files, fixture generators and full-data
reference drivers are preserved at
[fixed commit 2cb65ba5](https://github.com/hurry060215-tech/py-Scissor/tree/2cb65ba5af3e1ffac05226f66d0e0ddf1eb7dc50/reference), also named `validation-reference`.
This archive is separate from the current Python core and its source/wheel
packages. It contains no original LUAD matrices or individual real-data outputs.

Read the [historical audit and reproduction instructions](https://github.com/hurry060215-tech/py-Scissor/tree/2cb65ba5af3e1ffac05226f66d0e0ddf1eb7dc50/reference/README.md)
and [full-data oracle guide](https://github.com/hurry060215-tech/py-Scissor/tree/2cb65ba5af3e1ffac05226f66d0e0ddf1eb7dc50/reference/luad/README.md) there. R and a C++
compiler are needed only when deliberately regenerating independent development
references in that separate checkout. Ordinary installation, tests, examples and
all three Python models need neither.

The R-generated expected outputs are not rewritten during ordinary tests.
Missing fixtures or changed fixture hashes fail CI. Mathematical corrections and
upstream safety limits remain documented in [SEMANTICS](../docs/SEMANTICS.md).
