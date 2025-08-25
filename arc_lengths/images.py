from manim import *
import numpy as np
config.background_color = "#FEFAE1"   # any hex or "WHITE"/"BLACK"
config.renderer = "cairo"

class MajorMinorArcImage(Scene):

    def construct(self):
        # --- parameters ---
        R = 2.6
        start = 20*DEGREES          # where the minor arc starts
        delta = 120*DEGREES         # < pi so it's the "minor" arc

        def unit(a):
            return np.array([np.cos(a), np.sin(a), 0.0])

        # circle outline
        circle = Circle(radius=R, stroke_color=BLACK, stroke_width=6)

        # endpoints of the arcs
        A = R * unit(start)
        B = R * unit(start + delta)

        # minor (red) and major (green) arcs
        minor = Arc(radius=R, start_angle=start, angle=delta,
                    color="#F20000", stroke_width=14)
        major = Arc(radius=R, start_angle=start + delta, angle=TAU - delta,
                    color=GREEN, stroke_width=10).set_stroke(opacity=0.85)

        # chord and endpoint dots (helps show the same endpoints)
        chord = DashedLine(A, B, color=BLACK, stroke_width=3, dash_length=0.2)
        a_dot, b_dot = Dot(A, color=BLACK), Dot(B, color=BLACK)

        # labels near the midpoints of each arc
        minor_label = Text("minor arc", color=BLACK).scale(0.6)\
            .move_to(1.18 * R * unit(start + delta/2))
        major_label = Text("major arc", color=BLACK).scale(0.6)\
            .move_to(1.18 * R * unit(start + delta/2 + PI))

        # layer major first so minor sits on top
        self.add(circle, major, minor, chord, a_dot, b_dot, minor_label, major_label)
        self.wait(0.1)  # single frame only



class MajorMinorSectorImage(Scene):
    def construct(self):
        R = 2.6
        start = 20*DEGREES
        delta = 120*DEGREES   # < π → defines the minor sector

        def unit(a):
            return np.array([np.cos(a), np.sin(a), 0.0])

        # circle outline
        circle = Circle(radius=R, stroke_color=BLACK, stroke_width=6)

        # minor sector (red fill)
        minor_sector = Sector(
            radius=R,
            start_angle=start,
            angle=delta,
            fill_color="#F20000",
            fill_opacity=0.6,
            stroke_width=0
        )

        # major sector (green fill)
        major_sector = Sector(
            radius=R,
            start_angle=start + delta,
            angle=TAU - delta,
            fill_color=GREEN,
            fill_opacity=0.45,
            stroke_width=0
        )

        # radii boundaries
        OA = Line(ORIGIN, R*unit(start), color=BLACK, stroke_width=3)
        OB = Line(ORIGIN, R*unit(start+delta), color=BLACK, stroke_width=3)

        # dots at endpoints
        A = Dot(R*unit(start), color=BLACK)
        B = Dot(R*unit(start+delta), color=BLACK)

        # labels positioned roughly at centroid of each sector
        minor_label = Text("minor sector", color=BLACK).scale(0.6)\
            .move_to(R*0.55*unit(start + delta/2))
        major_label = Text("major sector", color=BLACK).scale(0.6)\
            .move_to(R*0.55*unit(start + delta/2 + PI))

        # layering: circle on top of fills, so the outline is visible
        self.add(major_sector, minor_sector, circle, OA, OB, A, B, minor_label, major_label)
        self.wait(0.1)   # single frame only