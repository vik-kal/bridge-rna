#!/bin/bash
#precompute.sh
source venv/bin/activate


python precompute/embed_osdr.py                              # NASA embeddings. Hours; resumable.
python precompute/build_projections.py  --skip-tsne                    # full-corpus projections. Hours, mostly 3-D t-SNE.
python precompute/fetch_archs4_meta.py                      # Earth metadata. ~35 s, needs network.
python precompute/validate_artifacts.py --mixing --quality  # gates the build


