# Third-party notices

Application source: Copyright 2026 TechTune Media Group, Apache License 2.0. See [LICENSE](LICENSE).

The application uses **Agrona 1.22.0**, licensed under Apache License 2.0. The native archive includes its original license from the [tagged upstream source](https://github.com/real-logic/agrona/tree/1.22.0).

The published Linux executable was compiled with **Oracle GraalVM 25.0.3+9.1**. Its embedded runtime components are governed by the [GraalVM Free Terms and Conditions](https://www.oracle.com/downloads/licenses/graal-free-license.html) and the applicable separate component terms. The archive preserves the compiler distribution's original `LICENSE.txt` and `legal` notices under `licenses/oracle-graalvm/`. This does not change the Apache-2.0 license on the application source. See the [GraalVM licensing FAQ](https://www.graalvm.org/faq/) for runtime-image distribution terms.

The executable dynamically links to the system's glibc and zlib; those shared libraries are not included in the archive.
