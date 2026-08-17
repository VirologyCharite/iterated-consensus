.PHONY: test lint upload nox wc clean

test:
	uv run pytest

lint:
	uv run ruff check .

# Requires PyPI credentials, e.g. UV_PUBLISH_TOKEN set in the environment
# (see https://docs.astral.sh/uv/guides/publish/). Runs lint and tests first.
upload: lint test clean
	uv build
	uv publish

nox:
	uv run noxfile.py

wc:
	@wc -l $$(find src tests -name '*.py' | sort)

clean:
	rm -rf dist build *.egg-info
	find src tests -name '__pycache__' -type d -exec rm -rf {} +
	rm -rf .pytest_cache .ruff_cache
