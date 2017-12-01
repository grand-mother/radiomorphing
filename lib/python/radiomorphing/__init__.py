"""A Python package for fast radio traces from air showers

Documentation coming sooon ...
"""

# Export the package public functions
from .scaling import scale
from .core import interpolate, process

# Register all modules
__all__ = ["core", "frame", "interpolation", "scaling", "utils"]
