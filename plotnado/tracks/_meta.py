"""Metaclass utilities for track model field shaping."""

from __future__ import annotations

import copy
from typing import Any, get_args, get_origin

from pydantic import BaseModel, Field
from pydantic_core import PydanticUndefined

_ModelMetaclass = type(BaseModel)


def _resolve_model_type(annotation: Any) -> type[BaseModel] | None:
    """Resolve a concrete Pydantic model class from an annotation."""
    if isinstance(annotation, type) and issubclass(annotation, BaseModel):
        return annotation

    origin = get_origin(annotation)
    if origin is None:
        return None

    args = get_args(annotation)
    if origin is list or origin is tuple or origin is dict:
        return None

    for arg in args:
        resolved = _resolve_model_type(arg)
        if resolved is not None:
            return resolved
    return None


class TrackMeta(_ModelMetaclass):
    """Inject aesthetics model fields into track subclasses for discoverability."""

    def __new__(
        mcs,
        name: str,
        bases: tuple[type, ...],
        namespace: dict[str, Any],
        **kwargs: Any,
    ) -> type:
        annotations = dict(namespace.get("__annotations__", {}))

        explicit_native_fields = {k for k in annotations if k != "aesthetics"}
        inherited_native_fields: set[str] = set()
        inherited_flattened_fields: set[str] = set()
        for base in bases:
            inherited_native_fields |= set(
                getattr(base, "__track_native_fields__", set())
            )
            inherited_flattened_fields |= set(
                getattr(base, "__track_flattened_aesthetics_fields__", set())
            )

        native_fields = explicit_native_fields | inherited_native_fields
        aesthetics_annotation = annotations.get("aesthetics")
        aesthetics_class = _resolve_model_type(aesthetics_annotation)
        flattened_fields: set[str] = set(inherited_flattened_fields)

        if aesthetics_class is not None:
            for field_name, field_info in aesthetics_class.model_fields.items():
                if field_name in native_fields:
                    continue

                annotations[field_name] = field_info.annotation

                if field_info.default is not PydanticUndefined:
                    namespace[field_name] = copy.deepcopy(field_info.default)
                elif field_info.default_factory is not None:
                    namespace[field_name] = Field(
                        default_factory=field_info.default_factory
                    )

                flattened_fields.add(field_name)

        namespace["__annotations__"] = annotations
        cls = super().__new__(mcs, name, bases, namespace, **kwargs)
        setattr(cls, "__track_native_fields__", frozenset(native_fields))
        setattr(
            cls,
            "__track_flattened_aesthetics_fields__",
            frozenset(flattened_fields),
        )
        setattr(cls, "__track_aesthetics_class__", aesthetics_class)
        return cls
