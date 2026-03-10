import coolbox.api as cb

def format_lim(lim):
    return int(lim) if float(lim) % 1 == 0 else f"{lim:.2f}"


def plot_text_range(self, ax, ymin, ymax, gr: cb.GenomeRange):
    ydelta = ymax - ymin

    # set min max
    ymax_print = format_lim(ymax)
    ymin_print = format_lim(ymin)
    small_x = 0.01 * gr.length

    if self.properties.get("data_range_location", 'left') == "right":
        ax.text(
            gr.end,
            ymax - ydelta * 0.2,
            f"[{ymin_print} - {ymax_print}]",
            horizontalalignment="right",
            verticalalignment="top",
        )
    else:
        # by default show the data range
        ax.text(
            gr.start,
            ymax - ydelta * 0.2,
            f"[{ymin_print} - {ymax_print}]",
            horizontalalignment="left",
            verticalalignment="top",
        )


def plot_label(self):
    if not self.properties.get("label_on_track"):
        if hasattr(self, "label_ax") and self.label_ax is not None:
            self.label_ax.text(
                0.15,
                0.5,
                self.properties["title"],
                horizontalalignment="left",
                size="large",
                verticalalignment="center",
            )

