from manim import *
from manim import CapStyleType, LineJointType, rate_functions as rf
from manim import Sector, ValueTracker
from manim import Text
import numpy as np


def make_box(size: tuple[float, float, float],
             center: tuple[float, float, float],
             color: str = "#c06") -> Cube:
    lx, ly, lz = size

    box = Cube(side_length=1.0)
    box.stretch_to_fit_depth(ly)
    box.stretch_to_fit_height(lz)
    box.stretch_to_fit_width(lx)

    # Solid fill, no individual edges
    box.set_fill(color=color, opacity=1.0)
    box.set_stroke(width=0)

    box.move_to(center)
    return box

class Compound(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        self.add(axes)

        self.set_camera_orientation(phi=75*DEGREES, theta=45*DEGREES)

        eps = 1e-3  # tiny offset

        box1 = make_box((4, 2, 2), (0, 0, 0))
        # lift the top block very slightly so contact faces don't coincide
        box2 = make_box((2, 2, 2), (-1, 0, 2 + eps))

        compound = VGroup(box1, box2)

        # Optional outline
        compound.set_stroke(color=BLACK, width=1)

        self.add(compound)
        self.move_camera(
            phi=90 * DEGREES,
            theta=0 * DEGREES,
            run_time=2
        )
        self.wait(1)
        self.wait(2)


class CompoundViews(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()

        # Build L-shaped solid
        eps = 1e-3
        base = make_box((4, 2, 2), (0, 0, 0))
        top  = make_box((2, 2, 2), (-1, 0, 2 + eps))
        solid = VGroup(base, top)

        self.add(axes, solid)

        # Start in a 3/4 perspective
        self.set_camera_orientation(phi=70 * DEGREES, theta=45 * DEGREES)
        self.wait(1)

        # -------- FRONT VIEW --------
        front_label = Text("Front view").scale(0.6).to_edge(DOWN)
        self.play(
            LaggedStart(
                FadeIn(front_label),
                run_time=0.5
            )
        )
        self.move_camera(
            phi=90 * DEGREES,
            theta=0 * DEGREES,
            run_time=2
        )
        self.wait(1)

        # -------- RIGHT SIDE VIEW --------
        side_label = Text("Right side view").scale(0.6).to_edge(DOWN)
        self.play(Transform(front_label, side_label))
        self.move_camera(
            phi=90 * DEGREES,
            theta=90 * DEGREES,
            run_time=2
        )
        self.wait(1)

        # -------- TOP VIEW --------
        top_label = Text("Top view").scale(0.6).to_edge(DOWN)
        self.play(Transform(front_label, top_label))
        self.move_camera(
            phi=0 * DEGREES,
            theta=0 * DEGREES,
            run_time=2
        )
        self.wait(1)
