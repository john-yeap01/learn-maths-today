from manim import *
from manim import CapStyleType, LineJointType, rate_functions as rf
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
        # Background circle
        arc_white = Arc(radius=1, start_angle=0, angle=TAU - 1e-4)
        arc_white.set_stroke("#478978", width=2)
        arc_white.set_cap_style(CapStyleType.ROUND)
        arc_white.joint_type = LineJointType.ROUND
        self.play(FadeIn(arc_white))
        self.wait(0.2)

        # --- Angle tracker + dynamic arc ---
        theta = ValueTracker(0.0)  # radians

        def make_arc():
            a = theta.get_value()
            return (
                Arc(radius=1, start_angle=0, angle=a)
                .set_stroke("#111827", width=5)
                .set_cap_style(CapStyleType.ROUND)
            )

        arc = always_redraw(make_arc)
        self.play(Create(arc))

        # --- Label (Text, no LaTeX) ---
        def deg_label():
            deg = theta.get_value() * 180 / PI
            return Text(f"{deg:0.1f}°", font_size=36).to_corner(UL).set_color("#111827")

        label = always_redraw(deg_label)
        self.add(label)

        # ------------------------------------------------------------------
        # Smooth keyframes: up → down → up, but with softer easing
        # ------------------------------------------------------------------
        self.play(
            theta.animate.set_value(220*DEGREES),
            run_time=1.8,
            rate_func=rf.ease_in_out_sine
        )
        self.play(
            theta.animate.set_value(130*DEGREES),
            run_time=1.5,
            rate_func=rf.ease_in_out_sine
        )
        self.play(
            theta.animate.set_value(300*DEGREES),
            run_time=2.0,
            rate_func=rf.ease_in_out_sine
        )
        self.wait(0.5)

        # ------------------------------------------------------------------
        # Arbitrary function: continuous smooth wobble
        # ------------------------------------------------------------------
        def theta_of_alpha(alpha: float) -> float:
            base = 180*DEGREES
            return base + (40*DEGREES) * np.sin(2*np.pi*alpha)

        self.play(
            UpdateFromAlphaFunc(theta, lambda m, a: m.set_value(theta_of_alpha(a))),
            run_time=3.0,
            rate_func=rf.linear  # linear in alpha, but sine makes it smooth
        )
        self.wait(0.5)

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


class Comparison(Scene):
    def construct(self):
        R = 1.0
        circle_fill   = "#517664"
        wedge_color   = "#F20000"
        outline_color = "#111827"
        radius_color  = "#8B95C9"
        text_color    = "#000000"
        EPS = 1e-3
        start = 0*DEGREES

        disk = Circle(radius=R).set_fill(circle_fill, 0.8).set_stroke(width=0)
        self.play(FadeIn(disk))
        self.wait(0.1)

        theta = ValueTracker(EPS)
        
        wedge = always_redraw(lambda:
            Sector(
                radius=R + EPS,             # tiny overlap to hide seams
                start_angle=start,
                angle=max(EPS, theta.get_value())
            ).set_fill(wedge_color, 1).set_stroke(width=0)
        )

        # # arc outline growing with the same theta
        # arc = always_redraw(lambda:
        #     Arc(
        #         radius=R,
        #         start_angle=start,
        #         angle=max(EPS, theta.get_value())
        #     ).set_stroke(outline_color, width=6).set_cap_style(CapStyleType.ROUND)
        # )

        # radius line that follows the tip of the wedge
        radius_line = always_redraw(lambda:
            Line(
                ORIGIN,
                R * np.array([
                    np.cos(start + max(EPS, theta.get_value())),
                    np.sin(start + max(EPS, theta.get_value())),
                    0.0
                ])
            ).set_stroke(radius_color, width=6).set_cap_style(CapStyleType.ROUND)
        )

        # angle readout (in degrees)
        readout = always_redraw(lambda:
            Text(
                f"{np.rad2deg(theta.get_value()):.0f}°",
            ).scale(0.7).to_corner(UL).set_fill(text_color, 1)
        )

        # add in correct z-order (wedge below radius/arc so outline is crisp)
        disk.z_index   = 0
        wedge.z_index  = 1
        # arc.z_index    = 2
        radius_line.z_index = 3

        self.add(wedge, radius_line, readout)
        # linear = constant speed; change to smootherstep if you want easing
        self.play(theta.animate.set_value(TAU - EPS), run_time=1.5, rate_func=linear)
        self.wait(0.3)

class AngleMarkerNumeric(Scene):
    def construct(self):
        R = 2.0
        start = 0*DEGREES
        EPS = 1e-3

        # main circle
        circle = Circle(radius=R).set_fill("#517664", 0.8).set_stroke(width=0)
        self.play(FadeIn(circle))

        theta = ValueTracker(EPS)

        # sweeping sector
        sector = always_redraw(lambda:
            Sector(
                radius=R + EPS,
                start_angle=start,
                angle=max(EPS, theta.get_value())
            ).set_fill("#F20000", 0.6).set_stroke(width=2)
        )

        # radius lines
        line1 = Line(ORIGIN, R*RIGHT).set_stroke("#111827", 3)   # fixed
        line2 = always_redraw(lambda:
            Line(
                ORIGIN,
                R * np.array([
                    np.cos(start + max(EPS, theta.get_value())),
                    np.sin(start + max(EPS, theta.get_value())),
                    0
                ])
            ).set_stroke("#111827", 3)
        )

        # central angle marker arc (small arc near origin)
        marker_r = 0.3
        angle_marker = always_redraw(lambda:
            Arc(
                radius=marker_r,
                start_angle=start,
                angle=max(EPS, theta.get_value())
            ).set_stroke("#000000", width=3)
        )

        # numeric label in degrees, positioned mid-arc
        label = always_redraw(lambda:
            Text(
                f"{np.rad2deg(theta.get_value()):.0f}°",
                color="#000000"
            ).scale(0.5).move_to(
                marker_r * np.array([
                    np.cos(start + theta.get_value()/2),
                    np.sin(start + theta.get_value()/2),
                    0
                ]) * 1.3   # push a bit outward
            )
        )



        self.add(sector, line1, line2, angle_marker, label)
        self.play(theta.animate.set_value(TAU - EPS), run_time=3, rate_func=linear)
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


class AngleSweepPulseTrace(Scene):
    def construct(self):
        # --- params / styling ---
        R = 1                      # circle radius
        marker_r = 0.15               # central angle marker radius
        start = 0*DEGREES
        EPS = 1e-3

        fill_col   = "#517664"
        sector_col = "#F20000"
        rim_col    = "#111827"
        marker_col = "#000000"
        label_col  = "#000000"

        # Pulse intensity (how big the temporary scale/width boost becomes)
        scale_amp      = 0.18   # ~+18% radius at the pulse peak
        width_base     = 5
        width_boost    = 10     # halo arc width at pulse peak
        halo_opacity   = 0.28

        # backdrop circle
        disk = Circle(radius=R).set_fill(fill_col, 0.8).set_stroke(width=0)
        self.play(FadeIn(disk))

        # driver for sweep
        theta = ValueTracker(EPS)

        # sweeping sector fill
        sector = always_redraw(lambda:
            Sector(
                radius=R + EPS,
                start_angle=start,
                angle=max(EPS, theta.get_value())
            ).set_fill(sector_col, 0.55).set_stroke(width=0)
        )

        # radius lines (fixed and moving)
        base_radius = Line(ORIGIN, R*RIGHT).set_stroke(rim_col, 3)
        unit = lambda a: np.array([np.cos(a), np.sin(a), 0.0])
        moving_radius = always_redraw(lambda:
            Line(ORIGIN, R * unit(start + max(EPS, theta.get_value())))
                .set_stroke(rim_col, 3)
        )

        # central angle marker + numeric label in degrees
        angle_marker = always_redraw(lambda:
            Arc(
                radius=marker_r,
                start_angle=start,
                angle=max(EPS, theta.get_value())
            ).set_stroke(marker_col, 3)
        )
        label = always_redraw(lambda:
            Text(f"{np.rad2deg(theta.get_value()):.0f}°", color=label_col)
            .scale(0.4)
            .move_to(1.8 * marker_r * unit(start + theta.get_value()/2))
        )

        # --- traced rim arc (normal) ---
        traced_arc = always_redraw(lambda:
            Arc(radius=R, start_angle=start, angle=max(EPS, theta.get_value()))
                .set_stroke(rim_col, width_base).set_cap_style(CapStyleType.ROUND)
        )

        # --- pulsing halo arc (inflates then returns to normal once during sweep) ---
        # bell-shaped pulse: 0 → peak at mid-sweep → 0
        def pulse(progress: float) -> float:
            # smooth bell (0..1): 4p(1-p); feel free to swap a different profile
            return 4 * progress * (1 - progress)

        halo_arc = always_redraw(lambda:
            (lambda prog:
                Arc(
                    radius=R * (1 + scale_amp * pulse(prog)),
                    start_angle=start,
                    angle=max(EPS, theta.get_value())
                )
                .set_stroke(rim_col,
                            width=int(np.interp(pulse(prog), [0,1], [width_base, width_boost])),
                            opacity=halo_opacity)
            )(theta.get_value() / TAU)
        )

        # tip dot that also pulses a bit
        tip_dot = always_redraw(lambda:
            (lambda prog:
                Dot(R * unit(start + theta.get_value()),
                    radius=np.interp(pulse(prog), [0,1], [0.05, 0.09]),
                    color=rim_col)
            )(theta.get_value() / TAU)
        )

        arc_length = always_redraw(lambda:
            Text(
                f"s = {R * theta.get_value():.2f}",   # s = R·θ (θ in radians)
                color="#111827"
            ).scale(0.6).to_corner(UL)                                   
        )

        # layering
        self.add(sector, base_radius, moving_radius,
                 traced_arc, halo_arc, tip_dot,
                 angle_marker, label, arc_length)

        self.play(theta.animate.set_value(TAU - EPS), run_time=3.0, rate_func=smootherstep)
        self.wait(0.7)

        # --- green glow effect at the end ---
        glow = Arc(radius=R, start_angle=start, angle=TAU)\
            .set_stroke(color=GREEN, width=15, opacity=0.6)\
            .set_cap_style(CapStyleType.ROUND)

        self.play(FadeIn(glow), run_time=0.3)
        self.wait(1.0)   # keep glowing
        self.play(FadeOut(glow), run_time=0.5)


        inner_glow = Arc(radius=marker_r, start_angle=start, angle=TAU)\
            .set_stroke(color=GREEN, width=15, opacity=0.6)\
            .set_cap_style(CapStyleType.ROUND)
        
        self.play(FadeIn(inner_glow), run_time=0.3)
        self.wait(1.0)   # keep glowing
        self.play(FadeOut(inner_glow), run_time=0.5)


        self.play(FadeOut(arc_length), run_time=0.2)

        arc_length_2 = Text()





class ShowArcs(Scene):
    def construct(self):
        # ------------------------------
        # Left: Arc of a parabola
        # ------------------------------
        parabola_func = lambda x: -0.6*(x**2) + 2

        # Full parabola (gray)
        parabola = FunctionGraph(parabola_func, x_range=[-2, 2], color=GRAY)

        # Highlighted arc (green)
        parabola_arc = FunctionGraph(parabola_func, x_range=[-1.7, 1.2], color=GREEN)

        # Endpoints
        A = Dot([-1.7, parabola_func(-1.7), 0])
        D = Dot([1.2, parabola_func(1.2), 0])

        label_A = Text("A").next_to(A, DL, buff=0.1).scale(0.5)
        label_D = Text("D").next_to(D, DR, buff=0.1).scale(0.5)

        title_left = Text("An arc of a parabola", color=GREEN).scale(0.5).next_to(parabola, UP)

        parabola_group = VGroup(parabola, parabola_arc, A, D, label_A, label_D, title_left).shift(LEFT*3)

        # ------------------------------
        # Right: Arc of a circle
        # ------------------------------
        circle = Circle(color=GRAY)

        # Highlighted arc (purple, top half)
        arc = Arc(radius=1, start_angle=PI, angle=PI, color=PURPLE)

        O = Dot([0,0,0])
        C = Dot(arc.point_from_proportion(0))
        B = Dot(arc.point_from_proportion(1))

        label_O = Text("O").next_to(O, DOWN, buff=0.1).scale(0.5)
        label_C = Text("C").next_to(C, DL, buff=0.1).scale(0.5)
        label_B = Text("B").next_to(B, DR, buff=0.1).scale(0.5)

        title_right = Text("An arc", color=PURPLE).scale(0.5).next_to(circle, UP)

        circle_group = VGroup(circle, arc, O, C, B, label_O, label_C, label_B, title_right).shift(RIGHT*3)

        # ------------------------------
        # Animation sequence
        # ------------------------------
        self.play(FadeIn(parabola_group))
        self.wait(1)
        self.play(FadeIn(circle_group))
        self.wait(2)