import coolbox.api as cb


def plot_text_range(self, ax, ymin, ymax, gr: cb.GenomeRange):
    ydelta = ymax - ymin

    # set min max
    format_lim = lambda lim: int(lim) if float(lim) % 1 == 0 else "{:.2f}".format(lim)
    ymax_print = format_lim(ymax)
    ymin_print = format_lim(ymin)
    small_x = 0.01 * gr.length

    if self.properties.get("data_range_location") == "right":
        ax.text(
            gr.end - small_x,
            ymax - ydelta * 0.2,
            "[ {} ~ {} ]".format(ymin_print, ymax_print),
            horizontalalignment="right",
            verticalalignment="top",
        )
    else:
        # by default show the data range
        ax.text(
            gr.start - small_x,
            ymax - ydelta * 0.2,
            "[ {} ~ {} ]".format(ymin_print, ymax_print),
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

