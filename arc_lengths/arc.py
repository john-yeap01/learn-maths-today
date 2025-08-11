from manim import *
from manim import CapStyleType, LineJointType

config.background_color = "#FEFAE1"   # any hex or "WHITE"/"BLACK"
config.renderer = "cairo"

class Hello(Scene):
    def construct(self):
        text = Text("Hello, Manim!")
        self.play(Write(text))
        self.wait()

class MakeSquare(Scene):
    def construct(self):
        square = Square().set_fill("#517664", 0.5)
        square.set_stroke("#111827", width=10)
        square.set_cap_style(CapStyleType.ROUND)   # extends ends a bit, closes the seam
        # keep sharp corners:
        square.joint_type = LineJointType.MITER
        self.play(Create(square))

class CreateCircle(Scene):
    def construct(self):
        circle = Circle()
        circle = Circle().set_fill(PINK, 0.5).set_stroke("#111827", width=3)
        self.play(Create(circle))

        square = Square().set_fill("#517664", 0.5).set_stroke("#111827", width=3)
        square.move_to(circle)
        self.play(Transform(circle, square))

        self.wait()

class Arcs(Scene):
    def construct(self):
        circle = Circle().set_fill("#517664", 0.7).set_stroke(width=0)
        self.play(FadeIn(circle))

        self.wait(0.5)

        outline = circle.copy().set_fill(opacity=0)
        outline.set_stroke("#111827", width=10)
        # outline.set_cap_style(CapStyleType.BUTT) # do not use BUTT which is the default
        outline.set_cap_style(CapStyleType.ROUND)
        outline.rotate(0.001)

        self.play(Create(outline), run_time=1.2)


