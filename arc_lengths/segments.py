from manim import *
from manim import CapStyleType, LineJointType, rate_functions as rf
from manim import Sector, ValueTracker
from manim import Text
import numpy as np

config.background_color = "#FEFAE1"   # any hex or "WHITE"/"BLACK"
config.renderer = "cairo"


class SegmentChordThenArc(Scene):
    def construct(self):
        # --- parameters ---
        R = 2.6
        start = 25*DEGREES        # first endpoint angle
        delta = 120*DEGREES       # < π to make the minor arc/segment

        def unit(a): return np.array([np.cos(a), np.sin(a), 0.0])

        # Geometry
        circle = Circle(radius=R, stroke_color=BLACK, stroke_width=6)
        A = R * unit(start)
        B = R * unit(start + delta)

        chord = Line(A, B, stroke_color="#111827", stroke_width=6)
        arc   = Arc(radius=R, start_angle=start, angle=delta,
                    color="#111827", stroke_width=10).set_stroke(opacity=0.9)

        # (Optional) a faint fill approximating the segment region for context
        # by converting the arc's polyline points + chord endpoints into a polygon
        arc_pts = [p for p in arc.get_points()]  # many small straight segments
        segment_fill = Polygon(*arc_pts, B, A) \
            .set_fill("#F20000", opacity=0.18).set_stroke(width=0)

        # Endpoints
        dotA, dotB = Dot(A, color=BLACK), Dot(B, color=BLACK)

        # Labels (will appear during highlights)
        chord_label = Text("line joining two points (chord)", font_size=36)\
            .next_to(chord.get_center(), DOWN).set_z_index(10).set_fill(BLACK)
        arc_label   = Text("arc between the same points", font_size=36)\
            .to_edge(RIGHT).set_fill(BLACK)

        # Intro
        title = Text("Circle segment = bounded by a chord and its arc",
                     font_size=42).to_edge(UP).set_fill(BLACK)
        self.play(Write(title))
        self.play(Create(circle), FadeIn(dotA), FadeIn(dotB))
        self.play(Create(chord), Create(arc))
        self.add(segment_fill)  # faint context (no animation to keep focus)

        # 1) Highlight the chord
        self.play(
            chord.animate.set_color("#F20000").set_stroke(width=12),
            rate_func=smooth, run_time=0.8
        )
        self.play(Write(chord_label))
        self.wait(0.4)

        # De‑emphasize chord
        self.play(chord.animate.set_color("#111827").set_stroke(width=6), FadeOut(chord_label), run_time=0.5)

        # 2) Highlight the arc
        self.play(
            arc.animate.set_color(GREEN).set_stroke(width=14, opacity=1.0),
            rate_func=smooth, run_time=0.8
        )
        self.play(FadeIn(arc_label, shift=DOWN))
        self.wait(0.6)

        # Settle to neutral colors and end on the idea
        self.play(
            arc.animate.set_color("#111827").set_stroke(width=10, opacity=0.9),
            FadeOut(arc_label),
            run_time=0.6
        )

        final_note = Text("A segment is the region between a chord and its arc.",
                          font_size=40).to_edge(DOWN).set_fill(BLACK)
        self.play(Write(final_note))
        self.wait(1.2)



class SegmentEquationMinimal(Scene):
    def construct(self):
        # --- geometry params ---
        R = 2.8
        start = 25*DEGREES
        delta = 120*DEGREES               # (< π) minor sector/segment

        def unit(a): return np.array([np.cos(a), np.sin(a), 0.0])

        # --- base lines in black ---
        circle = Circle(radius=R, stroke_color=BLACK, stroke_width=6).set_z_index(3)
        A = R*unit(start)
        B = R*unit(start + delta)

        OA = Line(ORIGIN, A, color=BLACK, stroke_width=4).set_z_index(4)
        OB = Line(ORIGIN, B, color=BLACK, stroke_width=4).set_z_index(4)
        chord = Line(A, B, color=BLACK, stroke_width=6).set_z_index(4)

        dotA, dotB = Dot(A, color=BLACK).set_z_index(5), Dot(B, color=BLACK).set_z_index(5)

        # --- fills for highlighting (start invisible) ---
        # sector (between OA, OB, and arc)
        sector_fill = Sector(radius=R, start_angle=start, angle=delta)\
            .set_fill(GREEN, 0).set_stroke(width=0).set_z_index(0)
        # triangle (between OA, OB, and chord)
        triangle_fill = Polygon(ORIGIN, A, B).set_fill(BLUE, 0).set_stroke(width=0).set_z_index(0)
        # segment (between chord AB and arc AB) — approximate with sampled arc polygon
        arc_pts = [R*unit(start + t) for t in np.linspace(0, delta, 64)]
        segment_fill = Polygon(*arc_pts).set_fill("#F20000", 0).set_stroke(width=0).set_z_index(0)
        # (Polygon closes from last point (B) back to first (A) -> chord)

        # --- lay down the drawing (all black lines first) ---
        self.play(Create(circle), Create(OA), Create(OB), Create(chord), FadeIn(dotA), FadeIn(dotB))

        # --- equation (split into parts so we can color words independently) ---
        t1 = Text("Area of segment", font_size=40, color=BLACK)
        t2 = Text("=", font_size=40, color=BLACK)
        t3 = Text("Area of sector", font_size=40, color=BLACK)
        t4 = Text("-", font_size=40, color=BLACK)  # use "−" if your font supports it
        t5 = Text("Area of triangle", font_size=40, color=BLACK)

        eq = VGroup(t1, t2, t3, t4, t5).arrange(RIGHT, buff=0.4).to_edge(DOWN)
        self.play(Write(eq))

        # --- highlight 1: SEGMENT ---
        self.add(segment_fill)                          # behind lines
        self.play(
            segment_fill.animate.set_fill("#F20000", 0.35),
            t1.animate.set_color("#F20000"),
            run_time=0.7
        )
        self.wait(0.9)
        self.play(
            segment_fill.animate.set_fill("#F20000", 0.0),
            t1.animate.set_color(BLACK),
            run_time=0.5
        )

        # --- highlight 2: SECTOR ---
        self.add(sector_fill)
        self.play(
            sector_fill.animate.set_fill(GREEN, 0.30),
            t3.animate.set_color(GREEN),
            run_time=0.9
        )
        self.wait(0.8)
        self.play(
            sector_fill.animate.set_fill(GREEN, 0.0),
            t3.animate.set_color(BLACK),
            run_time=0.5
        )

        # --- highlight 3: TRIANGLE ---
        self.add(triangle_fill)
        self.play(
            triangle_fill.animate.set_fill(BLUE, 0.30),
            t5.animate.set_color(BLUE),
            run_time=0.7
        )
        self.wait(0.8)
        self.play(
            triangle_fill.animate.set_fill(BLUE, 0.0),
            t5.animate.set_color(BLACK),
            run_time=0.5
        )

        self.wait(0.6)



class Triangle(Scene):
    def construct(self):
        # --- geometry params ---
        R = 2.8
        start = 25*DEGREES
        delta = 120*DEGREES               # (< π) minor sector/segment

        def unit(a): return np.array([np.cos(a), np.sin(a), 0.0])

        # --- base lines in black ---
        circle = Circle(radius=R, stroke_color=BLACK, stroke_width=6).set_z_index(3)
        A = R*unit(start)
        B = R*unit(start + delta)

        OA = Line(ORIGIN, A, color=BLACK, stroke_width=4).set_z_index(4)
        OB = Line(ORIGIN, B, color=BLACK, stroke_width=4).set_z_index(4)
        chord = Line(A, B, color=BLACK, stroke_width=6).set_z_index(4)

        dotA, dotB = Dot(A, color=BLACK).set_z_index(5), Dot(B, color=BLACK).set_z_index(5)

        # --- fills for highlighting (start invisible) ---
        # sector (between OA, OB, and arc)
        sector_fill = Sector(radius=R, start_angle=start, angle=delta)\
            .set_fill(GREEN, 0).set_stroke(width=0).set_z_index(0)
        # triangle (between OA, OB, and chord)
        triangle_fill = Polygon(ORIGIN, A, B).set_fill(BLUE, 0).set_stroke(width=0).set_z_index(0)
        # segment (between chord AB and arc AB) — approximate with sampled arc polygon
        arc_pts = [R*unit(start + t) for t in np.linspace(0, delta, 64)]
        segment_fill = Polygon(*arc_pts).set_fill("#F20000", 0).set_stroke(width=0).set_z_index(0)
        # (Polygon closes from last point (B) back to first (A) -> chord)


        angle_marker = Angle(OA,OB,radius=0.4).set_stroke(BLACK)

        # TRIANGLE 
        triangle_fill = Polygon(ORIGIN, A, B).set_fill(BLUE, 0).set_stroke(width=0).set_z_index(0)
        t5 = Text("Area of triangle = ½·r²·sin θ", font_size=40, color=BLACK).to_edge(DOWN)

        # --- lay down the drawing (all black lines first) ---
        self.play(Create(circle), Create(OA), Create(OB), Create(chord), FadeIn(dotA), FadeIn(dotB), FadeIn(angle_marker))


        self.wait(0.6)

        self.add(triangle_fill)
        self.play(
            triangle_fill.animate.set_fill(BLUE, 0.30),
            t5.animate.set_color(BLUE),
            run_time=0.7
        )
        self.wait(0.8)
        self.play(
            triangle_fill.animate.set_fill(BLUE, 0.0),
            t5.animate.set_color(BLACK),
            run_time=0.5
        )

        self.wait(0.6)




class SegmentNumericWorkedExample(Scene):
    def construct(self):
        # --- Parameters ---
        R = 2.0
        start_deg = 20
        end_deg   = 120
        start = start_deg * DEGREES
        theta_deg = end_deg - start_deg
        theta = theta_deg * DEGREES
        center = LEFT * 4  # left side of the screen

        # --- Base circle ---
        circle = Circle(radius=R).set_stroke(BLACK, 3).move_to(center)
        self.play(FadeIn(circle))

        # --- Sector (20° to 120°) ---
        sector = Sector(radius=R, start_angle=start, angle=theta)
        sector.move_arc_center_to(center)
        sector.set_fill("#FFD166", opacity=0.6).set_stroke("#FFD166", 0)
        self.play(FadeIn(sector))

        # --- Segment (arc + chord) ---
        arc = Arc(radius=R, start_angle=start, angle=theta).move_arc_center_to(center)
        A, B = arc.get_start(), arc.get_end()
        samples = 80
        arc_points = [arc.point_from_proportion(t/samples) for t in range(samples+1)]
        pts = arc_points + [arc_points[0]]
        segment = Polygon(*pts).set_fill("#06D6A0", 0.6).set_stroke("#06D6A0", 0)
        chord = Line(A, B).set_stroke(BLACK, 2)
        self.play(FadeIn(segment), Create(chord))

        # --- Radius markers ---
        OA = Line(center, A).set_stroke(BLACK, 2)
        OB = Line(center, B).set_stroke(BLACK, 2)
        self.play(Create(OA), Create(OB))

        # --- Angle marker ---
        angle_marker = Angle(OA, OB, radius=0.7, other_angle=False).set_stroke(BLACK)
        angle_label = Text("θ").scale(0.6).next_to(angle_marker, UP, buff=0.1).set_color(BLACK)
        self.play(Create(angle_marker), FadeIn(angle_label))

        # --- Equation text (left corner) ---
        eqn = Text("Area of segment = Area of sector - Area of triangle").set_color(BLACK)
        eqn.scale(0.9)
        eqn.to_corner(DOWN + LEFT)
        self.play(Write(eqn))

        # ===========================================================
        # --- Numeric worked example on the right side ---
        # ===========================================================
        example1 = Text("Given R = 2, θ = 100°").scale(0.8).set_color(BLACK)
        example1.to_edge(RIGHT, buff=4.0).shift(UP*3)
        self.play(Write(example1))
        self.wait(0.4)

        step1 = Text("Area of sector = 1/2 R² θ (in radians)").scale(0.8).set_color(BLACK)
        step1.next_to(example1, DOWN, aligned_edge=LEFT)
        self.play(Write(step1))
        self.wait(0.4)

        step2 = Text("= 1/2 × 2² × (100 × π/180)").scale(0.8).set_color(BLACK)
        step2.next_to(step1, DOWN, aligned_edge=LEFT)
        self.play(TransformMatchingShapes(step1.copy(), step2))
        self.wait(0.4)

        A_sector = 0.5 * R**2 * np.deg2rad(theta_deg)
        step3 = Text(f"= {A_sector:.3f}").scale(0.8).set_color(BLACK)
        step3.next_to(step2, DOWN, aligned_edge=LEFT)
        self.play(TransformMatchingShapes(step2.copy(), step3))
        self.wait(0.4)

        step4 = Text("Area of triangle = 1/2 R² sin θ").scale(0.8).set_color(BLACK)
        step4.next_to(step3, DOWN, aligned_edge=LEFT)
        self.play(Write(step4))
        self.wait(0.4)

        step5 = Text("= 1/2 × 2² × sin(100°)").scale(0.8).set_color(BLACK)
        step5.next_to(step4, DOWN, aligned_edge=LEFT)
        self.play(TransformMatchingShapes(step4.copy(), step5))
        self.wait(0.4)

        A_triangle = 0.5 * R**2 * np.sin(np.deg2rad(theta_deg))
        step6 = Text(f"= {A_triangle:.3f}").scale(0.8).set_color(BLACK)
        step6.next_to(step5, DOWN, aligned_edge=LEFT)
        self.play(TransformMatchingShapes(step5.copy(), step6))
        self.wait(0.4)

        A_segment = A_sector - A_triangle
        final = Text(f"Area of segment = {A_segment:.3f}").scale(0.7).set_color(RED)
        final.next_to(step6, DOWN, aligned_edge=LEFT)
        self.play(Write(final))
        self.wait(0.4)

        self.wait()