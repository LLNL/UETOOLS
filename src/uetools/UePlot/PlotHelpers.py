
class PlotHelpers:



    @staticmethod
    def _get_occupied_bboxes(ax, renderer, exclude=None):
        """
        Collect bbox-based obstacles that are *not* lines.
        Lines are handled separately with exact segment-vs-rectangle testing.
        """
        boxes = []
        exclude = set() if exclude is None else set(exclude)

        for artist in ax.texts:
            if artist in exclude or not artist.get_visible():
                continue
            try:
                bb = artist.get_window_extent(renderer)
                if bb.width > 0 and bb.height > 0:
                    boxes.append(bb)
            except Exception:
                pass

        for artist in ax.collections:
            if artist in exclude or not artist.get_visible():
                continue
            try:
                bb = artist.get_window_extent(renderer)
                if bb.width > 0 and bb.height > 0:
                    boxes.append(bb)
            except Exception:
                pass

        for artist in ax.patches:
            if artist in exclude or not artist.get_visible():
                continue
            try:
                bb = artist.get_window_extent(renderer)
                if bb.width > 0 and bb.height > 0:
                    boxes.append(bb)
            except Exception:
                pass

        return boxes

    @staticmethod
    def _all_line_clearance(ax, text_bbox, pad_px=2.0, exclude_lines=None):
        """
        Check the candidate text bbox against *all* visible lines in the axes.

        Returns
        -------
        intersects_any : bool
            True if the bbox intersects any line.
        min_dist : float
            Minimum display-space distance from the bbox to any line polyline.
        """
        from numpy import inf, asarray, isfinite, column_stack
        exclude_lines = set() if exclude_lines is None else set(exclude_lines)

        rx0 = text_bbox.x0 - pad_px
        ry0 = text_bbox.y0 - pad_px
        rx1 = text_bbox.x1 + pad_px
        ry1 = text_bbox.y1 + pad_px

        global_min_dist =inf

        for other_line in ax.lines:
            if other_line in exclude_lines or not other_line.get_visible():
                continue

            x = asarray(other_line.get_xdata(orig=False))
            y = asarray(other_line.get_ydata(orig=False))
            mask = isfinite(x) & isfinite(y)
            x = x[mask]
            y = y[mask]

            if len(x) < 2:
                continue

            pts = ax.transData.transform(column_stack([x, y]))

            # Hard reject if any segment intersects the label box
            for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
                if PlotHelpers._segment_intersects_rect(x0, y0, x1, y1, rx0, ry0, rx1, ry1):
                    return True, 0.0

            # Otherwise track nearest distance
            sample_points = [
                ((rx0 + rx1) * 0.5, (ry0 + ry1) * 0.5),
                (rx0, ry0), (rx0, ry1), (rx1, ry0), (rx1, ry1),
                ((rx0 + rx1) * 0.5, ry0), ((rx0 + rx1) * 0.5, ry1),
                (rx0, (ry0 + ry1) * 0.5), (rx1, (ry0 + ry1) * 0.5),
            ]

            for px, py in sample_points:
                for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
                    d = PlotHelpers._point_to_segment_distance(px, py, x0, y0, x1, y1)
                    if d < global_min_dist:
                        global_min_dist = d

        return False, global_min_dist

    @staticmethod
    def _bbox_overlap_area(bb1, bb2):
        x0 = max(bb1.x0, bb2.x0)
        y0 = max(bb1.y0, bb2.y0)
        x1 = min(bb1.x1, bb2.x1)
        y1 = min(bb1.y1, bb2.y1)
        if x1 <= x0 or y1 <= y0:
            return 0.0
        return (x1 - x0) * (y1 - y0)

    @staticmethod
    def _point_to_segment_distance(px, py, x0, y0, x1, y1):
        from numpy import hypot
        vx = x1 - x0
        vy = y1 - y0
        wx = px - x0
        wy = py - y0

        c1 = vx * wx + vy * wy
        if c1 <= 0:
            return hypot(px - x0, py - y0)

        c2 = vx * vx + vy * vy
        if c2 <= c1:
            return hypot(px - x1, py - y1)

        t = c1 / c2
        projx = x0 + t * vx
        projy = y0 + t * vy
        return hypot(px - projx, py - projy)


    @staticmethod
    def _segment_intersects_rect(x0, y0, x1, y1, rx0, ry0, rx1, ry1):
        # quick reject
        if max(x0, x1) < rx0 or min(x0, x1) > rx1 or max(y0, y1) < ry0 or min(y0, y1) > ry1:
            return False

        # endpoint inside
        if (rx0 <= x0 <= rx1 and ry0 <= y0 <= ry1) or (rx0 <= x1 <= rx1 and ry0 <= y1 <= ry1):
            return True

        def ccw(ax, ay, bx, by, cx, cy):
            return (cy - ay) * (bx - ax) > (by - ay) * (cx - ax)

        def seg_intersect(ax, ay, bx, by, cx, cy, dx, dy):
            return ccw(ax, ay, cx, cy, dx, dy) != ccw(bx, by, cx, cy, dx, dy) and \
                   ccw(ax, ay, bx, by, cx, cy) != ccw(ax, ay, bx, by, dx, dy)

        edges = [
            (rx0, ry0, rx1, ry0),
            (rx1, ry0, rx1, ry1),
            (rx1, ry1, rx0, ry1),
            (rx0, ry1, rx0, ry0),
        ]
        for ex0, ey0, ex1, ey1 in edges:
            if seg_intersect(x0, y0, x1, y1, ex0, ey0, ex1, ey1):
                return True
        return False


    @staticmethod
    def _line_bbox_clearance(ax, line, text_bbox, pad_px=2.0):
        """
        Returns:
            intersects: bool
            min_dist: minimum distance in display pixels from bbox boundary/center samples
                      to the line polyline
        """
        from numpy import isfinite, asarray, column_stack, inf
        x = asarray(line.get_xdata(orig=False))
        y = asarray(line.get_ydata(orig=False))
        mask = isfinite(x) & isfinite(y)
        x = x[mask]
        y = y[mask]

        if len(x) < 2:
            return False, inf

        pts = ax.transData.transform(column_stack([x, y]))

        rx0 = text_bbox.x0 - pad_px
        ry0 = text_bbox.y0 - pad_px
        rx1 = text_bbox.x1 + pad_px
        ry1 = text_bbox.y1 + pad_px

        # hard intersection test
        for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
            if PlotHelpers._segment_intersects_rect(x0, y0, x1, y1, rx0, ry0, rx1, ry1):
                return True, 0.0

        # distance from a few representative bbox points to the polyline
        sample_points = [
            ((rx0 + rx1) * 0.5, (ry0 + ry1) * 0.5),  # center
            (rx0, ry0), (rx0, ry1), (rx1, ry0), (rx1, ry1),  # corners
            ((rx0 + rx1) * 0.5, ry0), ((rx0 + rx1) * 0.5, ry1),
            (rx0, (ry0 + ry1) * 0.5), (rx1, (ry0 + ry1) * 0.5),
        ]

        min_dist = inf
        for px, py in sample_points:
            for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
                d = PlotHelpers._point_to_segment_distance(px, py, x0, y0, x1, y1)
                if d < min_dist:
                    min_dist = d

        return False, min_dist


    @staticmethod
    def _score_text_bbox(text_bbox, occupied_bboxes, ax_bbox, anchor_px, ax, line):
        """
        Lower is better.

        Hard constraints:
        - full text bbox must stay inside axes
        - must not intersect any existing text bbox
        - must not intersect any line in the axes
        """
        from numpy import inf, isfinite
        # Must stay fully inside plotting area
        if (
            text_bbox.x0 - 50 < ax_bbox.x0 or
            text_bbox.y0 - 50 < ax_bbox.y0 or
            text_bbox.x1 + 50 > ax_bbox.x1 or
            text_bbox.y1 + 50 > ax_bbox.y1
        ):
            return inf

        # Must not overlap existing labels / text / collections / patches
        score = 0.0
        for bb in occupied_bboxes:
            overlap = PlotHelpers._bbox_overlap_area(text_bbox, bb)
            if overlap > 0:
                return inf
            score += 1000.0 * overlap

        # Must not overlap *any* plotted line
        intersects_any, min_dist = PlotHelpers._all_line_clearance(ax, text_bbox, pad_px=2.0)
        if intersects_any:
            return inf

        # Prefer a bit more clearance from nearby lines
        if isfinite(min_dist):
            score += 1500.0 / max(min_dist, 1.0)

        # Mild penalty for wandering too far from the anchor point
        cx = 0.5 * (text_bbox.x0 + text_bbox.x1)
        cy = 0.5 * (text_bbox.y0 + text_bbox.y1)
        dx = cx - anchor_px[0]
        dy = cy - anchor_px[1]
        score += 0.03 * (dx * dx + dy * dy)

        return score

    @staticmethod
    def _local_normal_offsets_pts(ax, x, y, i):
        """
        Candidate offsets in points, biased perpendicular to the local curve direction.
        """
        from numpy import hypot
        n = len(x)
        if n < 2:
            return [(10, 8), (10, -8), (-10, 8), (-10, -8), (12, 0), (-12, 0)]

        i0 = max(0, i - 1)
        i1 = min(n - 1, i + 1)
        p = ax.transData.transform([[x[i0], y[i0]], [x[i1], y[i1]]])
        dx, dy = p[1] - p[0]
        norm = hypot(dx, dy)

        if norm == 0:
            return [(10, 8), (10, -8), (-10, 8), (-10, -8), (12, 0), (-12, 0)]

        tx, ty = dx / norm, dy / norm
        nx, ny = -ty, tx

        px_to_pt = 72.0 / ax.figure.dpi
        candidates = []

        for normal_mag in [10, 14, 18]:
            for tangential_mag in [-4, 0, 4]:
                for s in [1, -1]:
                    ox_px = s * normal_mag * nx + tangential_mag * tx
                    oy_px = s * normal_mag * ny + tangential_mag * ty
                    candidates.append((ox_px * px_to_pt, oy_px * px_to_pt))

        candidates.extend([(10, 8), (10, -8), (-10, 8), (-10, -8), (14, 0), (-14, 0)])
        return candidates

    @staticmethod
    def add_inline_label(ax, line, label, fontsize=None, color=None, bbox=True, fontweight=None):
        from numpy import asarray, isfinite, inf
        from matplotlib.pyplot import rcParams
        if label is None:
            return None

        x = asarray(line.get_xdata(orig=False))
        y = asarray(line.get_ydata(orig=False))
        mask = isfinite(x) & isfinite(y)
        x = x[mask]
        y = y[mask]

        if len(x) == 0:
            return None

        fig = ax.figure
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()

        occupied = PlotHelpers._get_occupied_bboxes(ax, renderer, exclude={line})
        ax_bbox = ax.get_window_extent(renderer)

        if fontweight is None:
            fontweight = "normal"

        # Avoid the exact endpoint; try a few points a bit upstream.
        idx_candidates = []
        for frac in [0.98, 0.95, 0.92, 0.88, 0.84, 0.78]:
            i = int(frac * (len(x) - 1))
            if 0 <= i < len(x):
                idx_candidates.append(i)
        idx_candidates = list(dict.fromkeys(idx_candidates))

        best = None
        best_score = inf

        for i in idx_candidates:
            anchor_px = ax.transData.transform((x[i], y[i]))
            offset_candidates = PlotHelpers._local_normal_offsets_pts(ax, x, y, i)

            for dx_pts, dy_pts in offset_candidates:
                ha = "left" if dx_pts >= 0 else "right"

                tmp = ax.annotate(
                    label,
                    xy=(x[i], y[i]),
                    xytext=(dx_pts, dy_pts),
                    textcoords="offset points",
                    ha=ha,
                    va="center",
                    rotation=0,
                    fontsize=fontsize if fontsize is not None else rcParams["font.size"],
                    fontweight=fontweight,
                    color=color if color is not None else line.get_color(),
                    bbox=(
                        dict(
                            boxstyle="round,pad=0.15",
                            facecolor="white",
                            edgecolor="none",
                            alpha=0.75,
                        ) if bbox else None
                    ),
                    clip_on=False,
                    zorder=line.get_zorder() + 1,
                )

                fig.canvas.draw()
                bb = tmp.get_window_extent(renderer)
                tmp.remove()

                score = PlotHelpers._score_text_bbox(bb, occupied, ax_bbox, anchor_px, ax=ax, line=line)

                if isfinite(score) and score < best_score:
                    best_score = score
                    best = dict(i=i, dx_pts=dx_pts, dy_pts=dy_pts, ha=ha)

        if best is None:
            return None

        return ax.annotate(
            label,
            xy=(x[best["i"]], y[best["i"]]),
            xytext=(best["dx_pts"], best["dy_pts"]),
            textcoords="offset points",
            ha=best["ha"],
            va="center",
            rotation=0,
            fontsize=fontsize if fontsize is not None else rcParams["font.size"],
            fontweight=fontweight,
            color=color if color is not None else line.get_color(),
            bbox=(
                dict(
                    boxstyle="round,pad=0.15",
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.75,
                ) if bbox else None
            ),
            clip_on=False,
            zorder=line.get_zorder() + 1,
        )



