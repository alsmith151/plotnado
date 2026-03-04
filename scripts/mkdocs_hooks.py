"""MkDocs hooks for build-time warning hygiene."""

from __future__ import annotations

import logging
import sys
import warnings


def on_startup(**_: object) -> None:
    warnings.filterwarnings(
        "ignore",
        message=r".*parseString.*deprecated.*",
    )
    warnings.filterwarnings(
        "ignore",
        message=r".*resetCache.*deprecated.*",
    )
    warnings.filterwarnings(
        "ignore",
        message=r".*enablePackrat.*deprecated.*",
    )

    class _SuppressExternalHtmlNoise(logging.Filter):
        def filter(self, record: logging.LogRecord) -> bool:
            message = record.getMessage()
            if " Div at " in message and "unclosed" in message and "closing implicitly" in message:
                return False
            return True

    logging.getLogger().addFilter(_SuppressExternalHtmlNoise())

    class _FilteredStderr:
        def __init__(self, wrapped):
            self._wrapped = wrapped

        def write(self, text: str) -> int:
            if (
                " Div at " in text
                and "unclosed" in text
                and "closing implicitly" in text
            ):
                return len(text)
            return self._wrapped.write(text)

        def flush(self) -> None:
            self._wrapped.flush()

        def isatty(self) -> bool:
            return self._wrapped.isatty()

        @property
        def encoding(self):
            return self._wrapped.encoding

    sys.stderr = _FilteredStderr(sys.stderr)
