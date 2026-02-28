# Testing python script to get familiar with Manim's 3D functions and quirks
# Scenes for the basic prisms
from manim import *
from manim import CapStyleType, LineJointType, rate_functions as rf
from manim import Sector, ValueTracker
from manim import Text
import numpy as np


# axes = ThreeDAxes()
# self.set_camera_orientation(phi=75*DEGREES, theta=45*DEGREES)
# self.begin_ambient_camera_rotation(rate=0.2)

class RectPrism(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        self.set_camera_orientation(phi=75*DEGREES, theta=45*DEGREES)
        self.begin_ambient_camera_rotation(rate=0.2)

        prism = Prism((3,2,4))

        self.add(axes, prism)

        self.begin_ambient_camera_rotation(rate=0.2)

        self.wait(4)

        self.stop_ambient_camera_rotation()


class CreateCylinder(ThreeDScene):
    def construct(self):

        axes = ThreeDAxes()
        self.set_camera_orientation(phi=75*DEGREES, theta=45*DEGREES)
        self.begin_ambient_camera_rotation(rate=0.2)

        cylinder = Cylinder(2, 5, [1, 2, 0])

        self.add(cylinder)
        self.wait(10)



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
        self.play(h.animate.set_value(3.0), run_time=4)
        self.wait(13)


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







# USING THE PRISM OBJECT TO EXTRUDE
class CylinderExtrude (ThreeDScene):

    def construct(self):
        axes = ThreeDAxes()
        self.set_camera_orientation(phi=75*DEGREES, theta=45*DEGREES)
        self.begin_ambient_camera_rotation(rate=0.6)

        side = 2.0
        h = ValueTracker(0.01)  # start almost flat

        cylinder = always_redraw(
            lambda: Cylinder(radius=1, height=h.get_value(), direction=[1,2,0] )
                    .set_fill(BLUE, 0.6).set_stroke(WHITE, 1)
                    .move_to([0, 0, h.get_value()/2])
        )

        self.add(axes, cylinder)
        self.play(h.animate.set_value(4.0), run_time=7)
        self.wait(1)



from manim import *

class TriangleAreaText(Scene):
    def construct(self):
        # Geometry values
        base_len = 4
        height_len = 3
        area_val = 0.5 * base_len * height_len

        # Triangle points
        A = LEFT * base_len / 2
        B = RIGHT * base_len / 2
        C = UP * height_len

        # Triangle
        triangle = Polygon(A, B, C, color=BLUE)
        triangle.set_fill(BLUE, opacity=0.5)

        # Base
        base = Line(A, B, color=YELLOW)
        base_label = Text(f"base = {base_len}", font_size=28).next_to(base, DOWN)

        # Height (perpendicular)
        foot = np.array([0, 0, 0])
        height = DashedLine(C, foot, color=RED)
        height_label = Text(f"height = {height_len}", font_size=28).next_to(height, RIGHT)

        # Area text (plain text)
        area_text = Text(
            f"Area = 1/2 × base × height = {area_val}",
            font_size=32
        ).to_edge(UP)

        # Animations
        self.play(Create(triangle))
        self.wait(0.3)

        self.play(Create(base), FadeIn(base_label))
        self.wait(0.3)

        self.play(Create(height), FadeIn(height_label))
        self.wait(0.3)

        self.play(FadeIn(area_text))
        self.wait(2)
