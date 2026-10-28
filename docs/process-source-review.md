# Process source review

The layout work used three kinds of input, with different qualifications.

**Original KIT reference.** The user-supplied `NY_CREATES` generator and design
report describe a **research-grade NY CREATES NbN process reference**. They
provide nominal geometry and the die/NbN/bump-pad/passivation stream map.
They do not establish a foundry-qualified minimum-width/spacing deck. The
bundled research profile preserves that distinction: nominal dimensions are
`reference_geometry`, while executable minimum checks are labeled
`engineering_assumption`. Default straight, turn and launch primitive regions
were compared directly against the original definition; see the
[zero-XOR record](../validation/layout-reference-primitives.json).

**Local SFQ material.** The SFQ project contains SeeQC customer rule PDFs and
an SFQ5ee KLayout layer setup. A viewer layer configuration is useful for
interpreting a stack but is not itself an executable fabrication rule deck.
Neither a multilayer SFQ stack nor an anodization-based junction process is
equivalent to the illustrative overlap-JJ stack shipped here.

**Local `spec2si-seeqc` repository.** Its process database and research DRC
work provide useful architectural lessons: retain per-record provenance,
separate drawn from fabricated area, distinguish rule transcription from
qualification, preserve unresolved stream datatypes, and describe contacts
explicitly. There are separate original and conservative Rev.12 rule sets;
conditional contacts, step restrictions and other interpretations cannot be
flattened into unconditional width/space tests. These lessons informed this
repository's process contract and validation boundaries.

The customer PDFs, their rule-table transcriptions, local PDK code and derived
SeeQC mask tables are **not bundled in this public repository**. The JJ example
has its own stream map and independently chosen values, all explicitly labeled
illustrative. It does not claim SeeQC compatibility or qualification. A local
SeeQC adapter can be developed separately against the reviewed rule-set
revision, once its stack/contact semantics and distribution conditions are
settled. Public literature alone should not override a selected local process
revision.

The next integration work is process-specific: qualified JJ primitives,
conditional DRC, bias/area calibration, inductance/capacitance extraction,
coupling/bias structures and extracted-circuit comparison with JoSIM. The
current repository supplies the process → geometry → independent checking
workflow and reproducible examples, rather than claiming those qualifications.
