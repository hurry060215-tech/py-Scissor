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
Only reviewed source should be executed; R reference code lives in a separate fixed development archive. Explicit Python
data-download helpers are development tools, never part of library runtime.
The project does not automatically download data, install packages or submit
analytics. OmicVerse and its integration tests belong to a separate contribution; it has
its own independent import behavior. It is not a dependency of this kernel.
