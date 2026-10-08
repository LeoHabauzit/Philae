import cadquery as cq
import numpy as np
from cadquery.vis import show
from microgen import Tpms
from microgen.shape import strut_lattice
from microgen.shape.surface_functions import gyroid


def linear(
    x,
    _,
    __,
):
    """Linearly graded offset."""
    min_offset = 0.6
    max_offset = 3.0

    print(f"max: {np.max(x)}, min: {np.min(x)}")
    # length = np.max(x) - np.min(x)
    return (min_offset + max_offset) / 2 + (max_offset - min_offset) / (
        np.max(x) - np.min(x)
    ) * x


def linear_graded_offset(
    x,
    _,
    __,
):
    """Piecewise continuous graded offset."""

    min_offset = 0.7698088443756911
    max_offset = 3.0
    x_start = 4.0
    x_end = 10.0
    x_abs = np.abs(x)

    return np.piecewise(
        x_abs,
        [
            x_abs <= x_start,
            (x_abs > x_start) & (x_abs < x_end),
            x_abs >= x_end,
        ],
        [
            min_offset,
            lambda x: (
                min_offset
                + (max_offset - min_offset) * (x - x_start) / (x_end - x_start)
            ),
            max_offset,
        ],
    )


geometry = Tpms(
    surface_function=gyroid,
    cell_size=2.0,
    offset=linear_graded_offset,
    repeat_cell=(18, 4, 4),
    resolution=30,
)

geometry.sheet.plot(color="white")
print(geometry.offset)
