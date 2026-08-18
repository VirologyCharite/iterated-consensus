#!/bin/sh
# Pre-commit hook: keep uv.lock in sync with pyproject.toml automatically.
#
# `uv lock` is idempotent -- a no-op if uv.lock already matches
# pyproject.toml (e.g. after a version bump). If it doesn't, this
# regenerates uv.lock, stages it (so it's included once you re-run the
# commit), and fails with a clear, specific message -- rather than letting
# a stale lockfile surface confusingly inside the pytest hook (see
# .pre-commit-config.yaml's comment on the pytest hook's --no-sync).
set -eu

before=$(git hash-object uv.lock)
uv lock -q
after=$(git hash-object uv.lock)

if [ "$before" != "$after" ]; then
    git add uv.lock
    echo "uv.lock was out of date (e.g. after a version bump in pyproject.toml)."
    echo "It has been regenerated and staged -- run 'git commit' again."
    exit 1
fi
