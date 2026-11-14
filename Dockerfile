# syntax=docker/dockerfile:1
ARG UV_IMAGE=ghcr.io/astral-sh/uv:python3.12-bookworm-slim

# ---------- builder: resolve + compile dependencies ----------
FROM ${UV_IMAGE} AS builder

ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PYTHON_DOWNLOADS=never \
    UV_PROJECT_ENVIRONMENT=/opt/venv

# Some deps (e.g. gtfreader via pyranges1) ship sdists only and need a C toolchain.
RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install third-party dependencies first (cached layer; only invalidated by lock changes).
# --frozen: install exactly what uv.lock pins without re-resolving.
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    --mount=type=bind,source=requirements.txt,target=requirements.txt \
    uv sync --frozen --no-dev --no-install-project --extra notebook --extra docs

# Install plotnado itself (non-editable so the venv is self-contained).
COPY pyproject.toml uv.lock requirements.txt README.md LICENSE MANIFEST.in ./
COPY plotnado ./plotnado
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-dev --no-editable --extra notebook --extra docs

# ---------- runtime ----------
FROM ${UV_IMAGE}

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH="/opt/venv/bin:$PATH" \
    HOME=/home/plotnado

# World-writable home so the image works under any `--user uid:gid`
# (Jupyter and matplotlib write config/runtime files to $HOME).
RUN mkdir -p /home/plotnado && chmod 1777 /home/plotnado

COPY --from=builder /opt/venv /opt/venv

WORKDIR /work

EXPOSE 8888

CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root", "--IdentityProvider.token="]
