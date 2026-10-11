# Safe reporting and research use

Do not place credentials, private patient data, individual clinical records or
access-controlled datasets in public issues, logs, notebooks or pull requests.
Prefer a minimal synthetic example for scientific or software bug reports.

For a confidential vulnerability, use GitHub's **Report a vulnerability** entry
under this repository's Security tab if it is available. If it is unavailable,
open only a non-sensitive request for a private reporting channel; do not post
exploit details or affected private data. No private email address or guaranteed
response time is advertised by this project.

The package is research pre-alpha and is not a clinical or diagnostic system.
Only reviewed source should be executed; R reference scripts and explicit data
download helpers are development tools, not part of ordinary library runtime.
The project does not automatically download data, install packages or submit
analytics. OmicVerse has its own independent import behavior; the included
integration tests use a network guard rather than changing that upstream code.
