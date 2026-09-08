# Experiment protocol

Each run writes `results/experiments/<experiment_id>/manifest.json` containing config, git commit, dataset manifest, input/output hashes, model version, random seed, software environment and UTC timestamp. Experiment identifiers must be immutable; an existing directory causes an explicit error.
