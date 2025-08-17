from manim import *

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

        arc = Arc(radius=1, start_angle=0, angle=PI)
        arc.set_stroke("#FFFFFF", 3)
        arc.move_arc_center_to(circle.get_center())
        self.play(Create(circle))
        self.play(Create(arc))

        # Transform into square::: 

        # square = Square().set_fill("#517664", 0.5).set_stroke("#111827", width=3)
        # square.move_to(circle)
        # self.play(Transform(circle, square))

        # self.wait()