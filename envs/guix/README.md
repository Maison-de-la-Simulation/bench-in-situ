

```bash
guix time-machine -C channels.scm -- shell -m manifest.scm
```

# Known issues

## Guix and python venv

When using python venv inside of a guix shell, there is a problem with certain pip installed libraries (e.g.: numpy, h5py).
These libraries fail with: `ImportError: libz.so.1: cannot open shared object file: No such file or directory`.

The solution is:
```bash
export LD_PRELOAD=$GUIX_ENVIRONMENT/lib/libz.so
```


