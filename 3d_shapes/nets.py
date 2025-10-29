from manim import *
import numpy as np
from typing import Tuple, Iterable, Dict
from numpy.typing import NDArray

def unit(v: Iterable[float]) :
    x = np.asarray(v, dtype=float).reshape(3)

    if x.size != 3:
        raise ValueError("Expected 3 components")
    
    n = float(np.linalg.norm(x))
    if n == 0.0:
        raise ValueError("zero length vector")
    
    return (x/n).astype(np.float64, copy=False)

# returns matrix cross form of a vector
def _skew(v: Iterable[float]):
    x = np.asarray(v, dtype=float)
    if x.size != 3:
        raise ValueError("Expected 3 components")
    x = x.reshape(3)
    # normalised: np.ndarray[Tuple[np.Any], np.dtype[np.float64]]= unit(x)
    x_, y_, z_ = map(float, x)
    skew = np.array([[0, -z_, y_], 
                                    [z_, 0, -x_],
                                    [-y_, x_, 0]], dtype=np.float64)
    return skew



# returns generalised rotation matrix given an axis and an angle around that axis
def rodrigues(axis: Iterable[float], theta: float):
    k = unit(axis)
    K = _skew(k)

    s = float(np.sin(theta))
    c = float(np.cos(theta))

    I = np.eye(3,dtype=np.float64)

    R = I + s*K + (1-c)*(K@K)

    return R



def make_square(side: float = 1.0, origin: Tuple[float, float]=(0.0, 0.0)):
    ox, oy = map(float, origin)
    s = float(side)
    V = np.array([
        [ox+0, oy+0, 0.0],
        [ox+s, oy+0, 0.0],
        [ox+s, oy+s, 0.0],
        [ox+0, oy+s, 0.0],
    ], dtype=np.float64)
    # edges as index pairs into V (counter-clockwise)
    E= np.array([[0,1],[1,2],[2,3],[3,0]], dtype=np.int32)
    return {"V": V, "E": E}


def se3_about_edge(p0: NDArray[np.float64],
                   p1: NDArray[np.float64],
                   theta: float):
    p0 = np.asarray(p0, dtype=np.float64).reshape(3)
    p1 = np.asarray(p1, dtype=np.float64).reshape(3)
    R = rodrigues(p1 - p0, theta)   # axis normalization handled inside rodrigues
    t = p0 - R @ p0                 # keep the entire line fixed
    return R, t


class HingePointDemo(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        axes = ThreeDAxes()
        self.add(axes)

        # hinge line (yellow): p0 -> p1 (x-axis segment)
        p0 = np.array([-1.0, -1.0, -1.0])
        p1 = np.array([1.0, 1.0, 1.0])
        hinge_line = Line3D(p0, p1, color=YELLOW)
        a = Dot3D(p0, color=YELLOW); b = Dot3D(p1, color=YELLOW)

        # point to rotate (not on the line)
        x0 = np.array([0.5, 0.6, 0.0])
        dot = Dot3D(x0, color=RED)

        theta = ValueTracker(0)

        def update_dot(mob: Mobject):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            mob.move_to(R @ x0 + t)
            return mob

        dot.add_updater(update_dot)
        self.add(hinge_line, a, b, dot)

        self.play(theta.animate.set_value(PI/2), run_time=2)  # 0→+90°
        self.play(theta.animate.set_value(0.0), run_time=2)   # back to 0
        dot.clear_updaters()

class HingeFaceStepA(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        axes = ThreeDAxes()
        self.add(axes)

        F = make_square(1.5, origin=(0.0, 0.0))  # 1.5×1.5 at (0,0)
        V, E = F["V"], F["E"]

        # 1) draw the face from V (CCW). Slight fill so edges are visible.
        face = Polygon(*[tuple(p) for p in V], color=BLUE, fill_opacity=0.2, stroke_width=2)
        self.add(face)

        # 2) highlight the hinge edge [0,1] in yellow (and mark endpoints)
        i0, i1 = map(int, E[0])          # first edge
        p0, p1 = V[i0], V[i1]
        hinge = Line3D(p0, p1, color=YELLOW)
        a = Dot3D(p0, color=YELLOW); b = Dot3D(p1, color=YELLOW)
        self.add(hinge, a, b)

        # 3) label vertices briefly to confirm ordering/orientation
        labels = VGroup(*[
            Text(str(i), font_size=24).move_to(p + np.array([0,0,0.02]))
            for i, p in enumerate(V)
        ])
        self.add(labels)

        self.wait(1)

class HingeFaceStepB(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        # original (flat) geometry
        F = make_square(1.5, origin=(0.0, 0.0))
        V0, E = F["V"].copy(), F["E"]          # keep an immutable copy V0
        i0, i1 = map(int, E[0])
        p0, p1 = V0[i0], V0[i1]                # hinge endpoints (world-fixed)

        # draw original face + hinge
        face = Polygon(*[tuple(p) for p in V0], color=BLUE, fill_opacity=0.25, stroke_width=2)
        hinge = Line3D(p0, p1, color=YELLOW)
        self.add(face, hinge, Dot3D(p0, color=YELLOW), Dot3D(p1, color=YELLOW))

        theta = ValueTracker(0.0)

        # updater: recompute rotated vertices and write back into the polygon
        def update_face(m: Polygon):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            V_rot = (R @ V0.T).T + t
            # Manim expects corners in order (no need to repeat first point for Polygon)
            m.set_points_as_corners([*map(tuple, V_rot)])
            return m

        face.add_updater(update_face)

        # animate: 0 → +90° → 0 (watch that the yellow edge stays exactly in place)
        self.play(theta.animate.set_value(PI/2), run_time=2)
        self.play(theta.animate.set_value(0.0), run_time=2)

        face.clear_updaters()
        self.wait(0.5)


class HingeFaceStepC(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        # original (flat) geometry
        F = make_square(1.5, origin=(0.0, 0.0))
        V0, E = F["V"].copy(), F["E"]          # keep an immutable copy V0
        i0, i1 = map(int, E[0])
        p0, p1 = V0[i0], V0[i1]                # hinge endpoints (world-fixed)

        # child geometry (we'll only show its hinge, no face)
        G = make_square(1.5, origin=(0.0, 0.0))
        V1, E1 = G["V"].copy(), G["E"].copy()
        j0, j1 = map(int, E1[2])               # child hinge indices 3->2
        j0, j1 = j1, j0                        # now direction is 3 -> 2
        q0, q1 = V1[j0], V1[j1]

        # (kept for minimal change; not strictly needed)
        u = unit(q1 - q0)
        k = unit(p1 - p0)

        # child hinge + dots (initial placement uses V1 directly)
        child_hinge = Line3D(V1[j0], V1[j1], color=GREEN)
        aG = Dot3D(V1[j0], color=GREEN)
        bG = Dot3D(V1[j1], color=GREEN)
        self.add(child_hinge, aG, bG)

        # draw original parent face + hinge
        face = Polygon(*[tuple(p) for p in V0], color=BLUE, fill_opacity=0.25, stroke_width=2)
        hinge = Line3D(p0, p1, color=YELLOW)
        self.add(face, hinge, Dot3D(p0, color=YELLOW), Dot3D(p1, color=YELLOW))

        theta = ValueTracker(0.0)

        # --- Updaters ---

        # parent folds about its hinge
        def update_face(m: Polygon):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            V_rot = (R @ V0.T).T + t
            m.set_points_as_corners([*map(tuple, V_rot)])
            return m
        face.add_updater(update_face)

        # child hinge rides along with the *same* parent transform
        def update_child_hinge(m: Line3D):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            P = (R @ V1[[j0, j1]].T).T + t      # uses broadcasting since the V1 points are rows 
            new_line = Line3D(P[0], P[1], color=GREEN)
            m.become(new_line)
            return m

        child_hinge.add_updater(update_child_hinge)

        def update_child_dot(m: Dot3D, idx: int):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            m.move_to((R @ V1[idx]) + t)
            return m
        aG.add_updater(lambda m: update_child_dot(m, j0))
        bG.add_updater(lambda m: update_child_dot(m, j1))

        # animate: 0 → +90° → 0
        self.play(theta.animate.set_value(PI/2), run_time=2)
        self.play(theta.animate.set_value(0.0), run_time=2)

        # cleanup
        face.clear_updaters()
        child_hinge.clear_updaters()
        aG.clear_updaters()
        bG.clear_updaters()
        self.wait(0.5)


class HingeFaceStepD(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70*DEGREES, theta=45*DEGREES)
        self.add(ThreeDAxes())

        # original (flat) geometry
        F = make_square(1.5, origin=(0.0, 0.0))
        V0, E = F["V"].copy(), F["E"]          # keep an immutable copy V0
        i0, i1 = map(int, E[0])
        p0, p1 = V0[i0], V0[i1]                # hinge endpoints (world-fixed)

        # child geometry (we'll only show its hinge, no face)
        G = make_square(1.5, origin=(0.0, 0.0))
        V1, E1 = G["V"].copy(), G["E"].copy()
        j0, j1 = map(int, E1[2])               # child hinge indices 3->2
        j0, j1 = j1, j0                        # now direction is 3 -> 2
        q0, q1 = V1[j0], V1[j1]

        # child face
        child_face = Polygon(*[tuple(p) for p in V1], color=GREEN, fill_opacity=0.25, stroke_width=2)

        # rotate the child face by -90 first (behind the scenes)
        R1, t1 = se3_about_edge(q0, q1, -PI/2)
        V_rot1 = (R1 @ V1.T).T + t1
        V1 = V_rot1
        child_face.set_points_as_corners([*map(tuple, V_rot1)])
        self.add(child_face)

        # (kept for minimal change; not strictly needed)
        u = unit(q1 - q0)
        k = unit(p1 - p0)

        # child hinge + dots (initial placement uses V1 directly)
        child_hinge = Line3D(V1[j0], V1[j1], color=GREEN)
        aG = Dot3D(V1[j0], color=GREEN)
        bG = Dot3D(V1[j1], color=GREEN)
        self.add(child_hinge, aG, bG)

        # draw original parent face + hinge
        face = Polygon(*[tuple(p) for p in V0], color=BLUE, fill_opacity=0.25, stroke_width=2)
        hinge = Line3D(p0, p1, color=YELLOW)
        self.add(face, hinge, Dot3D(p0, color=YELLOW), Dot3D(p1, color=YELLOW))


        # --- Rotation Updaters ---
        theta = ValueTracker(0.0)

        # parent folds about its hinge
        def update_face(m: Polygon):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            V_rot = (R @ V0.T).T + t
            m.set_points_as_corners([*map(tuple, V_rot)])
            return m
        face.add_updater(update_face)

        # child hinge rides along with the *same* parent transform
        def update_child_hinge(m: Line3D):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            P = (R @ V1[[j0, j1]].T).T + t      # uses broadcasting since the V1 points are rows 
            new_line = Line3D(P[0], P[1], color=GREEN)
            m.become(new_line)
            return m
        
        def update_child_face(m: Polygon):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            V_rot = (R @ V1.T).T + t
            m.set_points_as_corners([*map(tuple, V_rot)])
            return m

        child_hinge.add_updater(update_child_hinge)
        child_face.add_updater(update_child_face)

        def update_child_dot(m: Dot3D, idx: int):
            R, t = se3_about_edge(p0, p1, theta.get_value())
            m.move_to((R @ V1[idx]) + t)
            return m
        aG.add_updater(lambda m: update_child_dot(m, j0))
        bG.add_updater(lambda m: update_child_dot(m, j1))

        # animate: 0 → +90° → 0
        self.play(theta.animate.set_value(PI/2), run_time=2)
        self.play(theta.animate.set_value(0.0), run_time=2)

        # cleanup
        face.clear_updaters()
        child_hinge.clear_updaters()
        aG.clear_updaters()
        bG.clear_updaters()
        self.wait(0.5)