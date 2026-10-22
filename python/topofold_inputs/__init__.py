"""Optional input adapters for TopoFold (MDAnalysis required; not a dependency of the core ``topofold`` package)."""

from .ca_adapter import (
    ADAPTER_VERSION,
    AdapterError,
    CATrajectory,
    ResidueSpec,
    extract_ca,
    load_universe,
)

__all__ = ["ADAPTER_VERSION", "AdapterError", "CATrajectory", "ResidueSpec", "extract_ca", "load_universe"]
