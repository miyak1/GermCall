#!/usr/bin/env python3
"""Tiny YAML config reader used by the bash steps.

Usage:  cfg.py CONFIG.yaml KEY [DEFAULT]

KEY is a dotted path, e.g. "qc.min_dp".
  * scalar -> printed as-is (booleans as true/false)
  * list   -> one item per line
  * missing / null / empty -> DEFAULT (empty string if not given)
"""
import sys

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.stderr.write("ERROR: PyYAML is required (pip install pyyaml or conda install pyyaml)\n")
    sys.exit(2)


def main():
    if len(sys.argv) < 3:
        sys.stderr.write(__doc__)
        return 2
    path, key = sys.argv[1], sys.argv[2]
    default = sys.argv[3] if len(sys.argv) > 3 else ""

    with open(path) as fh:
        node = yaml.safe_load(fh) or {}

    for part in key.split("."):
        if isinstance(node, dict) and part in node:
            node = node[part]
        else:
            node = None
            break

    if node is None or node == "" or node == []:
        print(default)
    elif isinstance(node, bool):
        print("true" if node else "false")
    elif isinstance(node, (list, tuple)):
        for item in node:
            print(item)
    else:
        print(node)
    return 0


if __name__ == "__main__":
    sys.exit(main())
