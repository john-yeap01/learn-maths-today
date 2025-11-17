# Testing python script to get familiar with Manim's 3D functions and quirks
from manim import *
from manim import CapStyleType, LineJointType, rate_functions as rf
from manim import Sector, ValueTracker
from manim import Text
import numpy as np


# self.set_camera_orientation(phi=75*DEGREES, theta=45*DEGREES)


class Axes (ThreeDScene):

    def construct(self):
        axes = ThreeDAxes()
        self.set_camera_orientation(phi=75*DEGREES, theta=45*DEGREES)

        self.play(FadeIn(axes))
        # self.add(axes)


# 3d axes with a 2d face on the xy plane
class BaseOnAxes (ThreeDScene):

    def construct(self):
        axes = ThreeDAxes()
        self.set_camera_orientation(phi=75*DEGREES, theta=45*DEGREES)
        self.begin_ambient_camera_rotation(rate=0.2)

        square = Square(color= BLUE)

        self.play(FadeIn(axes, square))
        # self.add(axes)

        self.play(square.animate.scale(1.3))
        self.wait(2)
        self.play( square.animate.scale(0.5))
        self.wait(1)
        self.play(square.animate.scale(2))
        self.wait(2)

        self.stop_ambient_camera_rotation()


# USING THE PRISM OBJECT TO EXTRUDE
class PrismExtrude (ThreeDScene):

    def construct(self):
        axes = ThreeDAxes()
        self.set_camera_orientation(phi=75*DEGREES, theta=45*DEGREES)
        self.begin_ambient_camera_rotation(rate=0.2)

        side = 2.0
        h = ValueTracker(0.01)  # start almost flat

        prism = always_redraw(
            lambda: Prism(dimensions=[side, side, h.get_value()])
                    .set_fill(BLUE, 0.6).set_stroke(WHITE, 1)
                    .move_to([0, 0, h.get_value()/2])
        )

        self.add(axes, prism)
        self.play(h.animate.set_value(3.0), run_time=2)
        self.wait(1)


# class CubeExtrude (ThreeDScene):
     
#     def construct(self):
#         pass

        

# Test extrusion using vertices
class ExtrudeDemo(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.play(FadeIn(axes))

        base = Square().set_fill(BLUE, 0.6).set_stroke(WHITE, 1)
        top  = base.copy().set_fill(BLUE, 0.6).set_stroke(WHITE, 1)
        top.move_to(base)  # start co-planar

        # side faces will be built from corners; order matters (quad = 4 points)
        def sides():
            b = base.get_vertices()  # 4 points on z=0
            t = top.get_vertices()   # will move in z
            quads = [
                Polygon(b[0], b[1], t[1], t[0]),
                Polygon(b[1], b[2], t[2], t[1]),
                Polygon(b[2], b[3], t[3], t[2]),
                Polygon(b[3], b[0], t[0], t[3]),
            ]
            for q in quads:
                q.set_fill(BLUE, 0.35).set_stroke(WHITE, 1)
            return VGroup(*quads)

        side_group = always_redraw(sides)

        self.add(base, top, side_group)
        self.begin_ambient_camera_rotation(rate=0.15)

        # “Extrude”: lift the top face along +Z (OUT) while sides redraw
        self.play(top.animate.shift(OUT*2.5), run_time=2.5)

        self.wait(1)
        self.stop_ambient_camera_rotation()


