# Olla-DFT visual guide / Guía visual

Five original example images and a six-page PDF in each language. The final page shows local recovery validation from the benchmark. No new calculations were run and no numerical figure content was changed.

Each original PNG is copied byte for byte from the public software examples; `manifest.json` records SHA-256, captions, scope and versioned source links. The PDF is a presentation of those images, not a new scientific result. Existing figure labels may be in Spanish in both PDFs.

The examples use different inputs and methods. In particular, the LDA electronic band gap, scissor-corrected optical fit and phonon calculation are not interchangeable results. Check source conditions before reuse. Local recovery checks do not establish recovery after physical power loss or disk loss.

Author: Jorge Enrique González Sevilla. License: GPL-3.0-or-later, retained for these existing examples after the software transition to AGPL in 1.3.0. See ../../LICENSES/GPL-3.0.txt. Cite Olla-DFT, Quantum ESPRESSO and the pseudopotentials used.

Regenerate from the benchmark repository: `python tools/build_showcase.py --source ../olla-dft`.
