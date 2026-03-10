## Project Overview

PlotNado is a lightweight Python package for creating genome browser-style plots with a focus on simplicity, rich aesthetics, and performance. It provides a high-level API for composing complex genomic visualizations with multiple track types, advanced styling options, and automatic scaling and coloring features.

## Development Setup

- Use the .venv environment for development
- Activate with `source .venv/bin/activate` (Linux/Mac) or `.venv\Scripts\activate` (Windows)
- Install dependencies with `uv pip install ...`

## Code Style Guidelines

### General Principles
- Maintain clean, readable, and consistent code
- Keep methods small and focused on a single task
- Use type hints. Ensure they are accurate and comprehensive, including for function parameters and return types
- Ensure type hints are the most up to date i.e. list[str] instead of List[str] or str | None instead of Optional[str]
- Write unit tests for new functionality

### Comments
- Write comments in English only
- Add comments only when necessary to explain complex logic
- Don't comment obvious code
- Use docstrings for classes and public methods

### Naming Conventions
- Use descriptive names for classes, methods, and variables
- Follow Python naming conventions:
  - snake_case for methods and variables
  - CamelCase for classes
  - UPPER_CASE for constants