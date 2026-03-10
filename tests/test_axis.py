from plotnado.tracks import GenomicAxis, GenomicAxisAesthetics, GenomicRegion


class TestGenomicAxis:
    def test_axis_defaults_to_full_coordinates_without_chromosome(self, mock_ax):
        axis = GenomicAxis()
        region = GenomicRegion(chromosome="chr1", start=1_000_000, end=2_000_000)

        axis.plot(mock_ax, region)

        labels = [call.args[2] for call in mock_ax.text.call_args_list]
        assert labels
        assert all("M" not in str(label) for label in labels)
        assert all(str(label) != "chr1" for label in labels)
        assert all("," in str(label) and str(label).replace(",", "").isdigit() for label in labels)

    def test_axis_can_use_human_readable_labels(self, mock_ax):
        axis = GenomicAxis(aesthetics=GenomicAxisAesthetics(use_human_readable_labels=True))
        region = GenomicRegion(chromosome="chr1", start=1_000_000, end=2_000_000)

        axis.plot(mock_ax, region)

        labels = [str(call.args[2]) for call in mock_ax.text.call_args_list]
        assert any(label.endswith("M") for label in labels)