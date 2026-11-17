# Unreleased (after 2.2.5)

- Manual: how to open the macOS program (Gatekeeper warning for unsigned programs).
- The template `parameters` file has 100 cells and 1000 sites (it had 10 and 10).
- The macOS package has `libomp.dylib` writable, so `xattr -dr com.apple.quarantine .` works on the whole folder.
