# Phase 4.3.1 — Final encoder leakage hardening

Inspection shows that `VersionedEncoder` is configuration-defined and hashed from its immutable program/configuration specification. It has no `fit()` method and is not a trainable representation. The adversarial test therefore uses extreme future-only expression values and verifies that encoding those values cannot mutate the encoder hash or change the encoding of historical observations.

This establishes `frozen_encoder_future_mutation = UNIT_TESTED`. It does not establish leakage safety for a hypothetical trainable encoder; that status is `NOT_IMPLEMENTED`. The distinction is intentional and scientifically stronger than inventing an unimplemented fitting API.

The test suite remains at 57 collected/passing tests. Real-cohort forecasting is not implemented, Phase 5 remains blocked, simulation remains disabled, and no biological or causal conclusion follows from this software invariant.
