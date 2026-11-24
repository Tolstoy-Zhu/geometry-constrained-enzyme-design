# Local molecular renderer

- Package: 3dmol 2.5.5, the version used by py3Dmol 2.5.5.
- Source: https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js
- SHA-256: f7cc78921ae72e7623e89cdd111434f58c2efddd2ffda1cd212644b406fb8016
- License source: https://cdn.jsdelivr.net/npm/3dmol@2.5.5/LICENSE
- License SHA-256: 4c6eaaed856f3f28a3b1a98e74f4a8a71618de7d51ea4155c29f6f793bcef861
- Original license text retained in `3Dmol-LICENSE.txt` (BSD-3-Clause and incorporated component notices).

The local bundle is embedded before py3Dmol's guarded loader. The loader is pre-resolved, so the demo makes no CDN requests. No model weights are bundled.
