from manim import *
from manim import CapStyleType, LineJointType
from manim import Sector, ValueTracker

config.background_color = "#FEFAE1"   # any hex or "WHITE"/"BLACK"
config.renderer = "cairo"

class Hello(Scene):
    def construct(self):
        text = Text("Arcs")
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

class Arcs(Scene):
    def construct(self):
        # create circle
        circle = Circle().set_fill("#517664", 0.7).set_stroke(width=0)
        self.play(FadeIn(circle))
        self.wait(0.5)


        theta = PI

        # create ARC
        arc = Arc(radius=1, start_angle=0, angle=3.141)
        arc.set_stroke("#111827", width=3)
        arc.set_cap_style(CapStyleType.ROUND)
        arc.joint_type = LineJointType.ROUND
        self.play(Create(arc))

        # ROTATE THE WHOLE ARC
        # self.play(Rotate(arc, theta, about_point=circle.get_center()))



class Overlap(Scene):
    def construct(self):
        old = Circle().set_fill("#03091A", opacity=0.5)
        new = Circle().set_fill("#CEDD27", opacity=0.4)
        old.set_stroke(width=0)
        new.set_stroke(width=0)
         # optional safety:
        old.z_index = 1
        new.z_index = 0


        # order matters: new behind, old in front
        self.play(AnimationGroup(FadeIn(new), FadeIn(old)))

       

class Wipe(Scene):
    def construct(self):
        theta = ValueTracker(0)  # 90° so it’s visible

        cover = always_redraw(lambda:Sector(radius=1, start_angle=-90*DEGREES, 
                        angle=theta.get_value()).set_fill("#B71818", 1).set_stroke(width=0))
        self.add(cover)
        self.play(theta.animate.set_value(TAU), run_time=2, rate_func=linear)

        self.wait()


class Reveal(Scene):
    def construct(self):
        old = Circle().set_fill("#03091A", opacity=0.7)
        old.set_stroke(width=0)
        self.play(FadeIn(old)) # only the top color appears

        self.wait(0.2)

        # create ARC
        # arc = Arc(radius=1, start_angle=90*DEGREES, angle=359*DEGREES)
        # arc.set_stroke("#111827", width=3)
        # arc.set_cap_style(CapStyleType.ROUND)
        # arc.joint_type = LineJointType.ROUND
        # self.play(Create(arc))

        # USE NON ZERO epsilon
        EPS = 1e-3
        R = old.radius               # match circle size
        theta = ValueTracker(EPS)            # start slightly > 0 so it’s visible


        sector = always_redraw(lambda: 
            Sector(
                radius=R- EPS, 
                start_angle=90*DEGREES,     # top
                angle=max(EPS, theta.get_value())
            ).set_fill("#F802D3", 1).set_stroke(width=0)
        )

        arc = always_redraw(lambda:
            Arc(radius=1, start_angle=90*DEGREES, 
                angle=max(EPS, theta.get_value())
            ).set_stroke("#111827", width=3)
        )

        self.add(sector)
        self.play(theta.animate.set_value(TAU-EPS), run_time=1.2, rate_func=smootherstep)
        self.wait(0.2)

      

        