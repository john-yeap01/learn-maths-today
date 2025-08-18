from manim import *

class HandwrittenAngleArc(Scene):
    def construct(self):
        # --- Styling ---
        fs   = 54
        gap  = 0.35       # vertical gap around fraction bar

        title = Text("Angle–Arc Relationship", font_size=42).to_edge(UP)
        self.play(Write(title))
        self.wait(0.2)

        # --- Helper to build a stacked fraction (Text + Line) ---
        def make_fraction(numer_str, denom_str):
            numer = Text(numer_str, font_size=fs)
            denom = Text(denom_str, font_size=fs)

            width = max(numer.width, denom.width) + 0.3
            bar = Line(LEFT*width/2, RIGHT*width/2, stroke_width=6)

            group = VGroup(numer, bar, denom)
            numer.next_to(bar, UP, buff=gap)
            denom.next_to(bar, DOWN, buff=gap)
            return group, numer, bar, denom

        # --- Build both sides (initial form) ---
        lhs, lhs_num, lhs_bar, lhs_den = make_fraction("Angle (degrees)", "360°")
        rhs, rhs_num, rhs_bar, rhs_den = make_fraction("Arc length", "Circumference")
        eq  = Text("=", font_size=fs)

        equation = VGroup(lhs, eq, rhs).arrange(RIGHT, buff=0.8).scale(0.9)
        equation.move_to(ORIGIN)

        # --- Animate with default Write ---
        self.play(Create(lhs_bar))
        self.play(Write(lhs_num))
        self.play(Write(lhs_den))

        self.play(Write(eq))

        self.play(Create(rhs_bar))
        self.play(Write(rhs_num))
        self.play(Write(rhs_den))
        self.wait(0.5)

        # --- Rename: "Arc length" -> "s", "Circumference" -> "C" ---
        new_s = Text("s", font_size=fs).move_to(rhs_num)
        new_C = Text("C", font_size=fs).move_to(rhs_den)
        self.play(
            TransformMatchingShapes(rhs_num, new_s),
            TransformMatchingShapes(rhs_den, new_C),
        )
        # rebind handles to the new objects
        rhs_num = new_s
        rhs_den = new_C
        self.wait(0.3)

        # --- Rearrangement to isolate s:  s = (Angle/360°) × C ---
        # Targets: s on the left, equals, LHS fraction in the middle, "×", then C
        eq2   = Text("=", font_size=fs)
        times = Text("×", font_size=fs)

        # Build a temporary layout to grab final positions
        layout = VGroup(
            rhs_num.copy(),  # s
            eq2,
            lhs.copy(),      # the fraction Angle/360°
            times,
            rhs_den.copy()   # C
        ).arrange(RIGHT, buff=0.5).move_to(ORIGIN)

        s_target, eq2_target, lhs_target, times_target, C_target = layout

        # Prepare movements
        rhs_num.generate_target()
        rhs_num.target.move_to(s_target)

        lhs.generate_target()
        lhs.target.move_to(lhs_target)

        rhs_den.generate_target()
        rhs_den.target.move_to(C_target)

        # Fade out the old "=", replace with new one in position; remove rhs fraction bar
        self.play(
            FadeOut(eq),
            MoveToTarget(rhs_num),
            MoveToTarget(lhs),
            MoveToTarget(rhs_den),
            FadeOut(rhs_bar),
            run_time=1.2
        )
        # Draw "=" and "×" in their places
        eq2.move_to(eq2_target)
        times.move_to(times_target)
        self.play(Write(eq2), Write(times))

        self.wait(1)


from manim import *

class AreaSectorEquation(Scene):
    def construct(self):
        fs   = 54
        gap  = 0.35

        title = Text("Area Ratio for a Sector", font_size=42).to_edge(UP)
        self.play(Write(title))
        self.wait(0.2)

        # Helper: build a stacked fraction out of Text + Line
        def make_fraction(numer_str, denom_str):
            numer = Text(numer_str, font_size=fs)
            denom = Text(denom_str, font_size=fs)
            width = max(numer.width, denom.width) + 0.3
            bar = Line(LEFT*width/2, RIGHT*width/2, stroke_width=6)
            group = VGroup(numer, bar, denom)
            numer.next_to(bar, UP, buff=gap)
            denom.next_to(bar, DOWN, buff=gap)
            return group, numer, bar, denom

        # Build both sides (words first)
        lhs, lhs_num, lhs_bar, lhs_den = make_fraction("Angle (degrees)", "360°")
        rhs, rhs_num, rhs_bar, rhs_den = make_fraction("Sector area", "Area of circle")
        eq  = Text("=", font_size=fs)

        equation = VGroup(lhs, eq, rhs).arrange(RIGHT, buff=0.8).scale(0.9).move_to(ORIGIN)

        # Animate: default Write/Create
        self.play(Create(lhs_bar)); self.play(Write(lhs_num)); self.play(Write(lhs_den))
        self.play(Write(eq))
        self.play(Create(rhs_bar)); self.play(Write(rhs_num)); self.play(Write(rhs_den))
        self.wait(0.4)
        
        # Rename: "Sector area" -> "A_s", "Area of circle" -> "A_c"
        new_As = Text("A_s", font_size=fs).move_to(rhs_num)
        new_Ac = Text("A_c", font_size=fs).move_to(rhs_den)
        self.play(
            TransformMatchingShapes(rhs_num, new_As),
            TransformMatchingShapes(rhs_den, new_Ac),
        )
        rhs_num, rhs_den = new_As, new_Ac
        self.wait(0.2)

        # Rearrangement to isolate A_s:  A_s = (Angle/360°) × A_c
        eq2   = Text("=", font_size=fs)
        times = Text("×", font_size=fs)

        layout = VGroup(
            rhs_num.copy(),   # A_s
            eq2,
            lhs.copy(),       # Angle/360°
            times,
            rhs_den.copy()    # A_c
        ).arrange(RIGHT, buff=0.5).move_to(ORIGIN)

        As_target, eq2_target, lhs_target, times_target, Ac_target = layout

        rhs_num.generate_target(); rhs_num.target.move_to(As_target)
        lhs.generate_target();     lhs.target.move_to(lhs_target)
        rhs_den.generate_target(); rhs_den.target.move_to(Ac_target)

        self.play(
            FadeOut(eq),
            FadeOut(rhs_bar),       # drop old right-side fraction bar
            MoveToTarget(rhs_num),
            MoveToTarget(lhs),
            MoveToTarget(rhs_den),
            run_time=1.2
        )
        eq2.move_to(eq2_target); times.move_to(times_target)
        self.play(Write(eq2), Write(times))
        self.wait(1)
