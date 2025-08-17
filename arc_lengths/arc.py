from manim import *
from manim import CapStyleType, LineJointType
from manim import Sector, ValueTracker
from manim import Text
import numpy as np

config.background_color = "#FEFAE1"   # any hex or "WHITE"/"BLACK"
config.renderer = "cairo"
TICK_COLOR = "#E4572E"
TICK_W     = 3
TICK_LEN   = 0.18

class Hello(Scene):
    def construct(self):
        text = Text("Arcs", color="#478978")
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

class Arcs2(Scene):
    def construct(self):
        arc_white = Arc(radius=1, start_angle=0, angle=2*PI)
        arc_white.set_stroke("#478978", width=2)
        arc_white.set_cap_style(CapStyleType.ROUND)
        arc_white.joint_type = LineJointType.ROUND
        self.play(FadeIn(arc_white))
        self.wait(0.5)

        # create ARC
        arc = Arc(radius=1, start_angle=0, angle=PI)
        arc.set_stroke("#111827", width=5)
        arc.set_cap_style(CapStyleType.ROUND)
        arc.joint_type = LineJointType.ROUND
        self.play(Create(arc))

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

        # USE NON ZERO epsilon
        EPS = 1e-3
        R = old.radius               # match circle size
        theta = ValueTracker(EPS)            # start slightly > 0 so it’s visible


        sector = always_redraw(lambda: 
            Sector(
                radius=R- EPS, 
                start_angle=90*DEGREES,     # top
                angle=max(EPS, theta.get_value())
            ).set_fill("#478978", 1).set_stroke(width=0)
        )

        arc = always_redraw(lambda:
            Arc(radius=1, start_angle=90*DEGREES, 
                angle=max(EPS, theta.get_value())
            ).set_stroke("#8B95C9", width=3)
        )

        self.add(sector)
        self.play(theta.animate.set_value(TAU-EPS), run_time=1.2, rate_func=smootherstep)
        self.wait(0.2)

      



class ArcLinkedReveal(Scene):
    def construct(self):
        R = 1.0
        start = 90*DEGREES           # 12 o’clock
        EPS = 1e-3


        old = Circle(radius=R).set_fill("#8B95C9", 0.75).set_stroke(width=0)
        old.z_index = 1

        self.play(FadeIn(old))        # old on top

        theta = ValueTracker(EPS)     # drives BOTH the sector and the arc

        # the “reveal” wedge (paint it the same color as `new`)
        sector = always_redraw(lambda:
            Sector(
                radius=R + EPS,               # tiny overlap to avoid hairline seams
                start_angle=start,
                angle=max(EPS, theta.get_value())
            ).set_fill("#F20000", 1).set_stroke(width=0)
        )

        # the outline arc that grows with the same theta
        arc = always_redraw(lambda:
            Arc(
                radius=R,
                start_angle=start,
                angle=max(EPS, theta.get_value())
            ).set_stroke("#000000", width=5).set_cap_style(CapStyleType.ROUND)
        )

        text = always_redraw(lambda:
            Text(
                str(theta.get_value())
            ).to_corner(UL).set_fill("#000000", 1)                         
        )   

        # add both; they stay locked together as theta changes
        self.add(sector, arc, text)
        self.play(theta.animate.set_value(TAU - EPS), run_time=1.2, rate_func=smootherstep)

        
        self.wait(0.5)


class WhyRads(Scene):
    def construct(self):
        pass

class UnrollCircumference(Scene):
    def construct(self):
        # --- Parameters you can tweak ---
        R = 1
        fill_color = "#517664"   # circle fill
        line_color = "#111827"   # stroke + baseline
        stroke_w   = 6
        y_line     = -3.0        # vertical position of the unrolled line
        scale_factor = 0.4          #trim the unrolled squiggle a bit

        # --- Circle (disk) + the open circumference (Arc) ---
        disk = Circle(radius=R).set_fill(fill_color, 0.85).set_stroke(width=0)
        rim  = Arc(radius=R, start_angle=0, angle=TAU - 1e-4)  # open, so it won't double back
        rim.set_stroke(line_color, stroke_w).set_cap_style(CapStyleType.ROUND)

        self.play(FadeIn(disk))
        self.play(Create(rim), run_time=0.8)
        self.wait(0.2)

        # --- “Separate” the circumference (optional tiny pop-off) ---
        self.play(rim.animate.scale(1.05), run_time=0.8)

        # --- Unwrap mapping: (x,y) on the circle -> (s, 0) where s = R*theta ---
        length  = TAU * R + 1          # base line length
        start_x = -length / 2

        def unwrap(p: np.ndarray) -> np.ndarray:
            theta = np.arctan2(p[1], p[0])  # (-pi, pi]
            if theta < 0:
                theta += TAU                # [0, 2pi)
            s = R * theta * scale_factor                  # arc length from cut point (angle 0)
            return np.array([start_x + s, y_line, 0.0])

        # # Animate the unrolling
        # self.play(ApplyPointwiseFunction(unwrap, rim), run_time=1, rate_func=smoothererstep)

        # Snap to a perfect baseline (no brace/label/ticks; also no white ghost line)
        baseline = Line([start_x, y_line, 0], [start_x + length, y_line, 0])\
                    .set_stroke(line_color, stroke_w)\
                    .set_cap_style(CapStyleType.ROUND)

        # Replace instead of transform to avoid antialiasing artefacts
        # slow, smooth mapping onto the straight line:
        self.play(Transform(rim, baseline), run_time=1.2, rate_func=smooth)

        # cleanup to avoid the faint ghost seam:
        rim.set_stroke(line_color, stroke_w)          # reassert style
        self.wait(0.5)


# needs latex
class ArcLinkedReveal2(Scene):
    def construct(self):
        R = 1.0
        start = 90*DEGREES
        EPS = 1e-3

        old = Circle(radius=R).set_fill("#8B95C9", 0.75).set_stroke(width=0)
        old.z_index = 1
        self.play(FadeIn(old))

        theta = ValueTracker(EPS)

        # One object does both: fill (reveal) + stroke (outline)
        sector = always_redraw(lambda:
            Sector(
                radius=R + EPS,            # tiny overlap to avoid seam
                start_angle=start,
                angle=max(EPS, theta.get_value())
            )
            .set_fill("#F20000", 1)
            .set_stroke("#000000", width=3)
            .set_cap_style(CapStyleType.ROUND)
        )

        # Lightweight, live-updating numeric label (no re-creation each frame)
        deg = DecimalNumber(0, num_decimal_places=1)
        deg.add_updater(lambda m: m.set_value(np.degrees(theta.get_value())))
        deg.scale(0.6).next_to(old, DOWN)

        self.add(sector, deg)
        self.play(theta.animate.set_value(TAU - EPS), run_time=1.2, rate_func=linear)
        self.wait(0.2)