
class PlotHelpers:
    @staticmethod
    def _get_occupied_bboxes(ax, renderer, exclude=None):
        """Collect screen-space bounding boxes of visible artists."""
        boxes = []
        exclude = set() if exclude is None else set(exclude)

        for artist in ax.lines:
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

        for artist in ax.texts:
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
    def _bbox_overlap_area(bb1, bb2):
        x0 = max(bb1.x0, bb2.x0)
        y0 = max(bb1.y0, bb2.y0)
        x1 = min(bb1.x1, bb2.x1)
        y1 = min(bb1.y1, bb2.y1)
        if x1 <= x0 or y1 <= y0:
            return 0.0
        return (x1 - x0) * (y1 - y0)


    @staticmethod
    def _score_text_bbox(text_bbox, occupied_bboxes, ax_bbox, anchor_px):
        """
        Lower is better.
        Strongly penalize overlaps and going outside the axes.
        Mildly penalize distance from the anchor point.
        """
        score = 0.0

        for bb in occupied_bboxes:
            score += 1000.0 * PlotHelpers._bbox_overlap_area(text_bbox, bb)

        outside = 0.0
        if text_bbox.x0 < ax_bbox.x0:
            outside += ax_bbox.x0 - text_bbox.x0
        if text_bbox.y0 < ax_bbox.y0:
            outside += ax_bbox.y0 - text_bbox.y0
        if text_bbox.x1 > ax_bbox.x1:
            outside += text_bbox.x1 - ax_bbox.x1
        if text_bbox.y1 > ax_bbox.y1:
            outside += text_bbox.y1 - ax_bbox.y1
        score += 5000.0 * outside

        cx = 0.5 * (text_bbox.x0 + text_bbox.x1)
        cy = 0.5 * (text_bbox.y0 + text_bbox.y1)
        dx = cx - anchor_px[0]
        dy = cy - anchor_px[1]
        score += 0.2 * (dx * dx + dy * dy)

        return score


    @staticmethod
    def add_inline_label(ax, line, label, fontsize=None, color=None, bbox=False, fontweight='bold'):
        """
        Place a horizontal inline label near the end of a line,
        choosing an offset that minimizes overlap.
        """
        from numpy import asarray, isfinite, any, inf
        from matplotlib.pyplot import rcParams
        if label is None:
            return None

        x = asarray(line.get_xdata(orig=False))
        y = asarray(line.get_ydata(orig=False))

        mask = isfinite(x) & isfinite(y)
        if not any(mask):
            return None
        x = x[mask]
        y = y[mask]

        n = len(x)
        if n == 0:
            return None

        fig = ax.figure
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()

        occupied = PlotHelpers._get_occupied_bboxes(ax, renderer, exclude={line})
        ax_bbox = ax.get_window_extent(renderer)

        # Bias candidate anchor points toward the end of the line
        idx_candidates = []
        for di in [0, -1, -2, -3, -5, -8, -13]:
            i = n - 1 + di
            if 0 <= i < n:
                idx_candidates.append(i)
        idx_candidates = list(dict.fromkeys(idx_candidates))

        # Candidate offsets in points, biased to the right/end but with fallbacks
        offset_candidates = [
            (8, 0), (10, 6), (10, -6),
            (14, 10), (14, -10),
            (0, 8), (0, -8),
            (-8, 0), (-10, 6), (-10, -6),
            (-14, 10), (-14, -10),
            (18, 0), (-18, 0),
        ]

        best = None
        best_score = inf

        for i in idx_candidates:
            anchor_px = ax.transData.transform((x[i], y[i]))

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
                    color=color if color is not None else line.get_color(),
                    bbox=(
                        dict(
                            boxstyle="round,pad=0.15",
                            facecolor="white",
                            edgecolor="none",
                            alpha=0.75,
                        )
                        if bbox else None
                    ),
                    clip_on=True,
                    zorder=line.get_zorder() + 1,
                )

                fig.canvas.draw()
                bb = tmp.get_window_extent(renderer)
                score = PlotHelpers._score_text_bbox(bb, occupied, ax_bbox, anchor_px)

                if score < best_score:
                    best_score = score
                    best = dict(i=i, dx_pts=dx_pts, dy_pts=dy_pts, ha=ha)

                tmp.remove()

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
            color=color if color is not None else line.get_color(),
            fontweight=fontweight,
            bbox=(
                dict(
                    boxstyle="round,pad=0.15",
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.75,
                )
                if bbox else None
            ),
            clip_on=True,
            zorder=line.get_zorder() + 1,
        )

