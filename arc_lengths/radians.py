from manim import *
import numpy as np 

class RadianWrap(Scene):
    def construct(self):
        R = 2.0
        circle = Circle(radius=R, color="#111827").set_stroke(width=4)
        center = ORIGIN

        theta = ValueTracker(0.0)  # radians

        arc = always_redraw(lambda: Arc(
            radius=R, start_angle=0, angle=theta.get_value(),
            arc_center=center, color=YELLOW
        ).set_stroke(width=8))

        dot = always_redraw(lambda: Dot(
            point=circle.point_at_angle(theta.get_value()),
            color=YELLOW
        ))

        # Number line that labels with Text (no LaTeX)
        nl = NumberLine(
            x_range=[0, TAU, 1],
            length=6,
            include_numbers=True,
            label_constructor=Text,              # <- key change
            decimal_number_config={"num_decimal_places": 2},
        ).to_edge(DOWN)

        nl_dot = always_redraw(lambda: Dot(nl.n2p(theta.get_value()), color=YELLOW))
        connector = always_redraw(lambda: DashedLine(
            start=nl.n2p(theta.get_value()),
            end=dot.get_center(),
            dash_length=0.1
        ).set_stroke(opacity=0.6))

        # Labels with Text (use unicode π, θ)
        r_label = Text("r").next_to(Line(center, center + RIGHT*R), UP, buff=0.15)
        s_text  = always_redraw(lambda: Text(f"s = r·θ    (unit circle: s = θ)").scale(0.5).to_corner(UR))

        # π and 2π markers
        pi_lbl  = Text("π").scale(0.6).next_to(nl.n2p(PI), DOWN, buff=0.2)
        tau_lbl = Text("2π").scale(0.6).next_to(nl.n2p(TAU), DOWN, buff=0.2)

        # Small extra integer ticks (optional)
        ticks = VGroup()
        for k in range(1, 6):
            t = nl.get_tick(k, size=0.15)
            lbl = Text(str(k)).scale(0.35).next_to(nl.n2p(k), DOWN, buff=0.1)
            ticks.add(VGroup(t, lbl))

        self.play(Create(circle))
        self.play(FadeIn(r_label))
        self.play(Create(arc), FadeIn(dot))
        self.play(FadeIn(nl), FadeIn(nl_dot), FadeIn(connector), FadeIn(ticks), FadeIn(pi_lbl), FadeIn(tau_lbl), FadeIn(s_text))

        self.play(theta.animate.set_value(TAU), run_time=6, rate_func=linear)

        # Degree “conversion sticker” without LaTeX
        deg_sticker = VGroup(
            Rectangle(width=3.6, height=0.6, fill_opacity=0.9, color=RED, stroke_width=0),
            Text("× π/180").scale(0.7).set_color(WHITE)
        ).arrange().move_to(nl.get_center() + 2*UP + 2*RIGHT)

        self.play(FadeIn(deg_sticker, scale=0.5))
        self.wait(1)



class StackingRadii(Scene):
    def construct(self):
        R = 2.5
        circle = Circle(radius=R, color=BLUE_D).set_stroke(width=4)
        self.play(Create(circle), run_time=1.2)

        # radius line
        radius_line = Line(ORIGIN, ORIGIN + R * RIGHT, color=YELLOW).set_stroke(width=6)
        radius_label = Text("r = 1", color=YELLOW).scale(0.5)
        radius_label.next_to(radius_line, DOWN, buff=0.2)
        self.play(GrowFromCenter(radius_line), FadeIn(radius_label))

        # Title
        title = Text("Stacking radii along the rim → Radians").scale(0.5)
        title.to_edge(UP)
        self.play(FadeIn(title, shift=0.2*DOWN))

        # Value tracker for angle
        theta_tracker = ValueTracker(0.0)
        theta_readout = always_redraw(
            lambda: Text(f"θ = {theta_tracker.get_value():.3f} rad", color=TEAL_A).scale(0.5).to_corner(UR).shift(0.2*LEFT + 0.2*DOWN)
        )
        self.play(FadeIn(theta_readout))

        # Show 1 radian definition
        one_rad_arc = Arc(radius=0.6, start_angle=0, angle=1.0, arc_center=ORIGIN, color=RED).set_stroke(width=6)
        one_rad_label = Text("1 radian", color=WHITE).scale(0.5).next_to(one_rad_arc, RIGHT, buff=0.2)
        self.play(Create(one_rad_arc), FadeIn(one_rad_label, shift=0.2*UP))
        self.wait(0.5)
        self.play(FadeOut(one_rad_arc), FadeOut(one_rad_label))

        # Stack arcs
        start_angle = 0.0
        full_turn = TAU
        unit_angle = 1.0
        n_full = int(math.floor(full_turn / unit_angle))  # 6
        remainder = full_turn - n_full * unit_angle       # ~0.283

        def set_radius_to(angle):
            new_end = ORIGIN + R * np.array([math.cos(angle), math.sin(angle), 0])
            radius_line.become(Line(ORIGIN, new_end, color=YELLOW).set_stroke(width=6))

        for k in range(n_full):
            arc = Arc(R, start_angle, unit_angle, arc_center=ORIGIN, color=RED).set_stroke(width=8)
            self.play(
                Create(arc, run_time=0.7),
                UpdateFromAlphaFunc(radius_line, lambda m, a, a0=start_angle: set_radius_to(a0 + a*unit_angle), run_time=0.7),
                theta_tracker.animate.set_value(theta_tracker.get_value() + unit_angle),
            )
            label = Text(f"{k+1} rad").scale(0.4)
            end_angle = start_angle + unit_angle
            P = ORIGIN + (R+0.2) * np.array([math.cos(end_angle), math.sin(end_angle), 0])
            label.move_to(P)
            self.play(FadeIn(label), run_time=0.3)
            start_angle += unit_angle

        if remainder > 1e-6:
            rem_arc = Arc(R, start_angle, remainder, arc_center=ORIGIN, color=ORANGE).set_stroke(width=8)
            self.play(
                Create(rem_arc, run_time=0.6),
                UpdateFromAlphaFunc(radius_line, lambda m, a, a0=start_angle: set_radius_to(a0 + a*remainder), run_time=0.6),
                theta_tracker.animate.set_value(theta_tracker.get_value() + remainder),
            )

        # Finale
        finale_text = Text("Full circle = 2π rad ≈ 6.283").scale(0.6).to_edge(DOWN)
        self.play(circle.animate.set_stroke(TEAL_A, width=6))
        self.play(FadeIn(finale_text, shift=0.2*UP))
        self.wait(1.2)

        self.play(*[FadeOut(m) for m in self.mobjects])


from manim import *
import math
import numpy as np

class StackingRadiiWithAngle(Scene):
    def construct(self):
        R = 2.5
        circle = Circle(radius=R, color=BLUE_D).set_stroke(width=4)
        self.play(Create(circle), run_time=1.0)

        # --- Fixed start radius (reference line at angle = 0) ---
        start_radius = Line(ORIGIN, ORIGIN + R * RIGHT, color=WHITE).set_stroke(width=5)
        self.play(GrowFromCenter(start_radius), run_time=0.6)

        # --- Sweeping radius (moves as we add arcs) ---
        sweep_radius = Line(ORIGIN, ORIGIN + R * RIGHT, color=YELLOW).set_stroke(width=6)
        self.play(GrowFromCenter(sweep_radius), run_time=0.6)

        # Label r = 1 (conceptual unit circle; we just scale visually for clarity)
        r_label = Text("r = 1", color=YELLOW).scale(0.5)
        r_label.next_to(start_radius, DOWN, buff=0.2)
        self.play(FadeIn(r_label), run_time=0.3)

        # Title + live theta readout (in radians)
        title = Text("Stacking radii along the rim").scale(0.5)
        title.to_edge(UP)
        self.play(FadeIn(title, shift=0.2*DOWN), run_time=0.4)

        theta_tracker = ValueTracker(0.0)
        theta_readout = always_redraw(
            lambda: Text(f"θ = {theta_tracker.get_value():.3f} rad", color=TEAL_A)
                .scale(0.5)
                .to_corner(UR)
                .shift(0.2*LEFT + 0.2*DOWN)
        )
        self.play(FadeIn(theta_readout), run_time=0.3)

        # --- Central angle marker: filled wedge from 0 to current θ (small radius for visibility) ---
        # Using Sector for a neat filled angle marker
        angle_radius = 0.9
        angle_marker = always_redraw(
            lambda: Sector(
                arc_center=ORIGIN,
                radius=angle_radius,
                start_angle=0.0,
                angle=theta_tracker.get_value(),
                color=ORANGE,
                fill_opacity=0.5,
                stroke_width=0
            )
        )
        # A thin outline on top of the wedge
        angle_outline = always_redraw(
            lambda: Arc(
                radius=angle_radius,
                start_angle=0.0,
                angle=theta_tracker.get_value(),
                arc_center=ORIGIN,
                color=ORANGE
            ).set_stroke(width=3)
        )

        # Dynamic θ label near the middle of the wedge
        theta_arc_label = always_redraw(
            lambda: Text("angle θ", color=ORANGE).scale(0.45)
                .move_to(
                    ORIGIN + 1.15 * angle_radius
                    * np.array([math.cos(theta_tracker.get_value()/2),
                                math.sin(theta_tracker.get_value()/2), 0])
                )
        )

        self.add(angle_marker, angle_outline, theta_arc_label)

        # --- 1 radian definition quick show (optional) ---
        one_rad_arc = Arc(radius=0.6, start_angle=0, angle=1.0, arc_center=ORIGIN, color=RED).set_stroke(width=6)
        one_rad_label = Text("1 radian", color=WHITE).scale(0.5).next_to(one_rad_arc, RIGHT, buff=0.2)
        self.play(Create(one_rad_arc), FadeIn(one_rad_label, shift=0.1*UP), run_time=0.8)
        self.wait(0.3)
        self.play(FadeOut(one_rad_arc), FadeOut(one_rad_label), run_time=0.3)

        # --- Helper to rotate the sweeping radius to a target angle ---
        def set_radius_to(line, angle):
            new_end = ORIGIN + R * np.array([math.cos(angle), math.sin(angle), 0])
            line.become(Line(ORIGIN, new_end, color=line.get_color()).set_stroke(width=line.get_stroke_width()))

        # --- Tick marks at arc ends (purely visual) ---
        def tick_at_end(end_angle, length=0.18, color=WHITE):
            P = ORIGIN + R * np.array([math.cos(end_angle), math.sin(end_angle), 0])
            tangent = np.array([-math.sin(end_angle), math.cos(end_angle), 0])
            return Line(P - (length/2)*tangent, P + (length/2)*tangent, color=color).set_stroke(width=3)

        # --- Stack arcs of length 1 rad until full circle (2π) ---
        start_angle = 0.0
        full_turn = TAU         # 2π
        unit_angle = 1.0
        n_full = int(math.floor(full_turn / unit_angle))  # 6
        remainder = full_turn - n_full * unit_angle       # ~0.283185...

        arcs = VGroup()
        ticks = VGroup()

        for k in range(n_full):
            arc = Arc(
                radius=R,
                start_angle=start_angle,
                angle=unit_angle,
                arc_center=ORIGIN,
                color=RED
            ).set_stroke(width=8, opacity=0.95)
            arcs.add(arc)

            self.play(
                Create(arc, run_time=0.7),
                UpdateFromAlphaFunc(
                    sweep_radius,
                    lambda m, a, a0=start_angle: set_radius_to(m, a0 + a*unit_angle),
                    run_time=0.7
                ),
                theta_tracker.animate.set_value(theta_tracker.get_value() + unit_angle),
            )

            # Tick & tiny label at each whole radian
            end_angle = start_angle + unit_angle
            tick = tick_at_end(end_angle)
            ticks.add(tick)
            tick_label_pos = ORIGIN + (R + 0.25) * np.array([math.cos(end_angle), math.sin(end_angle), 0])
            tick_label = Text(f"{k+1} rad").scale(0.4).move_to(tick_label_pos)
            self.play(FadeIn(tick), FadeIn(tick_label), run_time=0.25)
            self.play(FadeOut(tick_label), run_time=0.15)

            start_angle += unit_angle

        if remainder > 1e-6:
            rem_arc = Arc(
                radius=R,
                start_angle=start_angle,
                angle=remainder,
                arc_center=ORIGIN,
                color=ORANGE
            ).set_stroke(width=8, opacity=0.95)
            arcs.add(rem_arc)

            self.play(
                Create(rem_arc, run_time=0.6),
                UpdateFromAlphaFunc(
                    sweep_radius,
                    lambda m, a, a0=start_angle: set_radius_to(m, a0 + a*remainder),
                    run_time=0.6
                ),
                theta_tracker.animate.set_value(theta_tracker.get_value() + remainder),
            )

        # --- Finale: recolor & print the takeaway ---
        finale_text = Text("Full circle = 2π rad ≈ 6.283").scale(0.6).to_edge(DOWN)
        self.play(circle.animate.set_stroke(TEAL_A, width=6))
        self.play(LaggedStart(*[arc.animate.set_color(TEAL_A) for arc in arcs], lag_ratio=0.05, run_time=0.8))
        self.play(FadeIn(finale_text, shift=0.2*UP), run_time=0.6)
        self.wait(1.0)

